###############
## GAM+MQR ####
###############
#' Multiple Quantile Regression Using Generalised Additive Models and Linear Quantile Regression
#'
#' This function fits multiple conditional linear quantile regression models to the residuals of
#' a generalised additive model using \code{mgcv} and with facilities for cross-validation.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory
#' variables. May optionally contain a collumn called "kfold" with
#' numbered/labeled folds and "Test" for test data. If \code{data} contains a column called
#' "gam_pred" then gam modelling will be skipped and this will be used in quantile regression. 
#' @param formala A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right passed to \code{gam()}
#' or \code{bam()} from \code{mgcv}.
#' @param formula_qr Formula for linear quantile regression model for GAM residuals. Term \code{gam_pred}
#' is the prediction from the above GAM may be included in this formula. If null, the "terms" of the GAM model
#' are used as features in linear quantile regression.
#' @param model_res2 If \code{TRUE} also model squared residuals of GAM using a GAM. Defaults to \code{FALSE}.
#' @param formula_res2 Formula for GAM to predict squared residuals.
#' @param quantiles The quantiles to fit models for.
#' @param gbm_params List of parameters to be passed to \code{fit.gbm()}.
#' @param cv_folds Control for cross-validation if not supplied in \code{data}.
#' @param Sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param ... Additional agruments passter to \code{gam()} or (\code{bam()}).
#' @param use_bam If \code{TRUE} (default) then GAM is fit using (\code{bam()}) in stead of \code{gam()}. \code{bam} is better suited to large datasets but not all \code{gam} model options are available with \code{bam}. See \code{bam()} documentation for further details.
#' @param w Weights on the contribution of data to model fit. See \code{gam()}.
#' @param SortLimits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @details The returned predictive quantiles and GAM predictions are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' @return Returns a \code{list} containing predictive quantiles (in a \code{MultiQR} object), GAM models, and deterministic predictions
#' from GAMs.
#' @keywords Quantile Regression
#' @export
qreg_gam <- function(data,
                     formula,
                     formula_qr=NULL,
                     model_res2=F,
                     formula_res2 = formula,
                     quantiles=c(0.25,0.5,0.75),
                     cv_folds=NULL,
                     use_bam=T,
                     # w=rep(1,nrow(data)), # Not working...
                     sort=T,sort_limits=NULL,
                     ...){
  
  ### Set-up Cross-validation
  TEST<-F # Flag for Training (with CV) AND Test output
  if("kfold" %in% colnames(data)){
    if(!is.null(cv_folds)){warning("Using column \"kfold\" from data. Argument \"cv_folds\" is not used.")}
    
    if("Test" %in% data$kfold){
      TEST<-T
      nkfold <- length(unique(data$kfold))-1
    }else{
      nkfold <- length(unique(data$kfold))
    }
  }else if(is.null(cv_folds)){
    data$kfold <- rep(1,nrow(data))
    nkfold <- 1
  }else{
    data$kfold <- sort(rep(1:cv_folds,length.out=nrow(data)))
    nkfold <- cv_folds
  }
  
  if(!"BadData" %in% colnames(data)){
    data[,BadData:=F]
  }
  
  # Special names...
  if("gam_resid" %in% names(data)){warning("data$gam_resid will be overwritten!")}
  if("gam_pred" %in% names(data)){warning("data$gam_pred will be overwritten!")}
  
  ## GAM for conditional expectation (and possible squared residuals)
  
  ## ADD FUNCTIONALITY TO USE PRE-FIT GAM
  
  OUTPUT_MODEL <- list(call=list(formula=formula,
                                 formula_qr=formula_qr,
                                 formula_res2 = if(model_res2){formula_res2}else{NULL}),
                       gam_pred=data.table::copy(data[,.(BadData,kfold,y=get(as.character(formula[[2]])))]),
                       gams = list(),
                       rqs = list())
  
  formula_res2 <- reformulate(attr(terms(formula_res2),"term.labels"),response = "gam_res2")
  for(fold in unique(data$kfold)){
    
    print(paste0("GAM, kfold=",fold))
    
    ## Fit gams using either gam() or bam()
    gam_fit_method <- if(use_bam){bam}else{gam}
    
    # GAM for conditional expectation
    OUTPUT_MODEL$gams[[fold]] <- gam_fit_method(data=data[kfold!=fold & kfold!="Test" & BadData==F,],
                                                formula = formula,...)
    # GAM for squared residuals
    if(model_res2){ 
      OUTPUT_MODEL$gams[[paste0(fold,"_r")]] <- gam_fit_method(data=cbind(data[kfold!=fold & kfold!="Test" & BadData==F,],data.table(gam_res2=residuals(OUTPUT_MODEL$gams[[fold]])^2)),
                                                               formula = formula_res2,...)
    }
    
    
    ## Out-of-sample cross-validation predictions
    OUTPUT_MODEL$gam_pred[kfold==fold,
                          gam_pred:=predict(OUTPUT_MODEL$gams[[fold]],newdata = data[kfold==fold,])]  
    
  }
  
  ## Out-of-sample cross-validation residuals
  OUTPUT_MODEL$gam_pred[,gam_resid:=y-gam_pred]
  data[,gam_resid := OUTPUT_MODEL$gam_pred$gam_resid]
  data[,gam_pred := OUTPUT_MODEL$gam_pred$gam_pred]
  
  ## Container for predictive quantiles
  predqs <- data.table(matrix(as.numeric(NA),ncol = length(quantiles), nrow = nrow(data)))
  colnames(predqs) <- paste0("q",100*quantiles)
  
  
  ## Qunatile regression
  if(!is.null(formula_qr)){
    ## Use user-specified model equation for quantile regression
    formula_qr <- reformulate(attr(terms(formula_qr),"term.labels"),response = "gam_resid")
    for(fold in unique(data$kfold)){
      print(paste0("MQR, kfold=",fold))
      for(i in 1:length(quantiles)){
        
        OUTPUT_MODEL$rqs[[fold]][[paste0("q",100*quantiles[i])]] <- 
          rq(formula_qr,
             tau = quantiles[i],
             data = data[kfold!=fold & kfold!="Test" & BadData==F,],
             method = "br")
        
        predqs[data$kfold==fold,i] <- OUTPUT_MODEL$gam_pred[kfold==fold,gam_pred] +
          predict.rq(OUTPUT_MODEL$rqs[[fold]][[paste0("q",100*quantiles[i])]],data[kfold==fold,])
      }
    }
  }else{
    ## QR with features from GAM
    for(fold in unique(data$kfold)){
      print(paste0("MQR, kfold=",fold))
      ## Get training Data
      train <- predict(OUTPUT_MODEL$gams[[fold]],newdata = data[kfold!=fold & kfold!="Test" & BadData==F,],type = "terms")
      if(model_res2){
        train2 <- predict(OUTPUT_MODEL$gams[[paste0(fold,"_r")]],newdata = data[kfold!=fold & kfold!="Test" & BadData==F,],type = "terms")
        # Only need to retrun smooth terms as linear terms included already...
        train2 <- train2[,grep("\\(",colnames(train2)),drop=F]
        colnames(train2) <- paste0(colnames(train2),"_r")
        train <- cbind(train,train2); rm(train2)
      }
      train <- cbind(as.data.table(train),data[kfold!=fold & kfold!="Test" & BadData==F,.(gam_resid)])
      
      ## Out-of-sample data
      test_cv <- as.data.table(predict(OUTPUT_MODEL$gams[[fold]],newdata = data[kfold==fold,],type = "terms"))
      if(model_res2){
        test_cv2 <- as.data.table(predict(OUTPUT_MODEL$gams[[paste0(fold,"_r")]],newdata = data[kfold==fold,],type = "terms"))
        test_cv2 <- test_cv2[,grep("\\(",colnames(test_cv2)),drop=F]
        colnames(test_cv2) <- paste0(colnames(test_cv2),"_r")
        test_cv <- cbind(test_cv,test_cv2); rm(test_cv2)
      }
      
      for(i in 1:length(quantiles)){
        ## Fit QR model
        OUTPUT_MODEL$rqs[[fold]][[paste0("q",100*quantiles[i])]] <- 
          rq(formula = gam_resid ~ .,
             tau = quantiles[i],
             data = train,
             method = "br")
        
        ## Make predictions
        predqs[data$kfold==fold,i] <- OUTPUT_MODEL$gam_pred[kfold==fold,gam_pred] +
          predict.rq(OUTPUT_MODEL$rqs[[fold]][[paste0("q",100*quantiles[i])]],
                     newdata=test_cv)
        
      }
    }
  }
  
  
  data[,gam_resid:=NULL]
  data[,gam_pred:=NULL]
  
  class(predqs) <- c("MultiQR",class(predqs))
  
  ### ADD DEFAULT MODEL FOR PREDICT FUNCTION
  OUTPUT_MODEL[["default_model"]] <- if("Test"%in%unique(data$kfold)){"Test"}else{unique(data$kfold)[1]}
  class(OUTPUT_MODEL) <- c("qreg_gam",class(OUTPUT_MODEL))
  
  if(sort){
    predqs <- SortQuantiles(data = predqs,Limits = sort_limits)
  }
  
  
  return(list(mqr_pred=predqs,
              fit_mods=OUTPUT_MODEL))
  
}


#' Predict from model based on Generalised Additive Model and Linear Quantile Regression
#'
#' This function predicts from multiple conditional linear quantile regression models of the residuals of
#' a generalised additive model fit using \code{mqr_qreg_gam}.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @export
predict.qreg_gam <- function(object,
                             newdata = NULL,
                             quantiles = NULL,
                             cv_fold=NULL, ### USE DEFAULT
                             sort = T,
                             sort_limits = NULL){
  
  ## Checks
  if(class(object)[1]!="qreg_gam"){stop("object of wrong class, expecting \"qreg_gam\"")}
  # Correct columns in newdata
  # Availability of quantile models
  
  ## Use default modeul unless specified
  if(is.null(cv_fold)){
    cv_fold <- object$default_model
  }
  
  if(is.null(quantiles)){
    quantiles <- as.numeric(gsub(pattern = "q",replacement = "",names(object$rqs[[cv_fold]])))/100
  }
  if(!all(quantiles %in% (as.numeric(gsub(pattern = "q",replacement = "",names(object$rqs[[cv_fold]])))/100))){
    stop("Models not availalbe for all requested quantiles.")
  }
  
  ## Initialise containers
  predqs <- data.table(matrix(as.numeric(NA),ncol = length(quantiles), nrow = nrow(newdata)))
  colnames(predqs) <- paste0("q",100*quantiles)
  
  OUTPUT <- list(gam_pred=predict(object$gams[[cv_fold]],newdata = newdata),
                 mqr_pred=NULL)
  newdata[,gam_pred:=OUTPUT$gam_pred]
  
  ## Qunatile regression
  if(!is.null(object$call$formula_qr)){
    ## Use user-specified model equation for quantile regression
    for(i in 1:length(quantiles)){
      predqs[,i] <- OUTPUT$gam_pred +
        predict.rq(object = object$rqs[[cv_fold]][[paste0("q",100*quantiles[i])]],
                   newdata = newdata)
    }
  }else{
    ## QR with features from GAM
    newdata_terms <- as.data.table(predict(object$gams[[cv_fold]],newdata = newdata,type = "terms"))
    if(!is.null(object$formula_res2)){
      newdata_terms2 <- as.data.table(predict(object$gams[[paste0(cv_fold,"_r")]],newdata = newdata,type = "terms"))
      newdata_terms2 <- newdata_terms2[,grep("\\(",colnames(newdata_terms2)),drop=F]
      colnames(newdata_terms2) <- paste0(colnames(newdata_terms2),"_r")
      newdata_terms <- cbind(newdata_terms,newdata_terms2); rm(newdata_terms2)
    }
    
    for(i in 1:length(quantiles)){
      ## Make predictions
      predqs[,i] <- OUTPUT$gam_pred +
        predict.rq(object = object$rqs[[cv_fold]][[paste0("q",100*quantiles[i])]],
                   newdata = newdata_terms)
      
    }
  }
  
  class(predqs) <- c("MultiQR",class(predqs))
  
  if(sort){
    predqs <- SortQuantiles(data = predqs,Limits = sort_limits)
  }
  
  
  OUTPUT$mqr_pred <- predqs
  
  return(OUTPUT)
  
} 











#' DEPRECIATED: Multiple Quantile Regression Using Generalised Additive Models and Linear Quantile Regression
#'
#' OLD VERSION SUPERCEDED BY \code{qreg_gam}. This function fits multiple conditional linear quantile regression models to the residuals of
#' a generalised additive model using \code{mgcv} and with facilities for cross-validation.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory
#' variables. May optionally contain a collumn called "kfold" with
#' numbered/labeled folds and "Test" for test data. If \code{data} contains a column called
#' "gam_pred" then gam modelling will be skipped and this will be used in quantile regression. 
#' @param formala A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right passed to \code{gam()}
#' or \code{bam()} from \code{mgcv}.
#' @param formula_qr Formula for linear quantile regression model for GAM residuals. Term \code{gam_pred}
#' is the prediction from the above GAM may be included in this formula.
#' @param model_res2 If \code{TRUE} also model squared residuals of GAM using a GAM. Defaults to \code{FALSE}.
#' @param formula_res2 Formula for GAM to predict squared residuals.
#' @param quantiles The quantiles to fit models for.
#' @param gbm_params List of parameters to be passed to \code{fit.gbm()}.
#' @param cv_folds Control for cross-validation if not supplied in \code{data}.
#' @param Sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param ... Additional agruments passter to \code{gam()} or (\code{bam()}).
#' @param use_bam If \code{TRUE} (default) then GAM is fit using (\code{bam()}) in stead of \code{gam()}. \code{bam} is better suited to large datasets but not all \code{gam} model options are available with \code{bam}. See \code{bam()} documentation for further details.
#' @param w Weights on the contribution of data to model fit. See \code{gam()}.
#' @param SortLimits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @details The returned predictive quantiles and GAM predictions are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' @return Returns a \code{list} containing predictive quantiles (in a \code{MultiQR} object), GAM models, and deterministic predictions
#' from GAMs.
#' @keywords Quantile Regression
#' @export
MQR_qreg_gam <- function(data,
                         formula,
                         formula_qr=NULL,
                         model_res2=F,
                         formula_res2 = formula,
                         quantiles=c(0.25,0.5,0.75),
                         CVfolds=NULL,
                         ...,
                         use_bam=T,
                         w=rep(1,nrow(data)),
                         Sort=T,SortLimits=NULL){
  
  ### Set-up Cross-validation
  TEST<-F # Flag for Training (with CV) AND Test output
  if("kfold" %in% colnames(data)){
    if(!is.null(CVfolds)){warning("Using column \"kfold\" from data. Argument \"CVfolds\" is not used.")}
    
    if("Test" %in% data$kfold){
      TEST<-T
      nkfold <- length(unique(data$kfold))-1
    }else{
      nkfold <- length(unique(data$kfold))
    }
  }else if(is.null(CVfolds)){
    data$kfold <- rep(1,nrow(data))
    nkfold <- 1
  }else{
    data$kfold <- sort(rep(1:CVfolds,length.out=nrow(data)))
    nkfold <- CVfolds
  }
  
  if(!"BadData" %in% colnames(data)){
    data[,BadData:=F]
  }
  
  
  ## GAM for conditional expectation (and possible squared residuals)
  
  ## Skip and only do quantile regression if gam_pred is provided
  if(is.null(data$gam_pred)){
    
    data[,gam_pred:=as.numeric(NA)]
    gams <- list()
    formula_res2 <- reformulate(attr(terms(formula_res2),"term.labels"),response = "gam_res2")
    for(fold in unique(data$kfold)){
      ## Fit gam
      print(paste0("GAM, kfold=",fold))
      
      if(use_bam){
        gams[[fold]] <- bam(data=data[kfold!=fold & kfold!="Test" & BadData==F,],
                            formula = formula,
                            weights = w,...)
        if(model_res2){
          gams[[paste0(fold,"_r")]] <- bam(data=cbind(data[kfold!=fold & kfold!="Test" & BadData==F,],data.table(gam_res2=residuals(gams[[fold]])^2)),
                                           formula = formula_res2,
                                           weights = w,...)
        }
      }else{
        gams[[fold]] <- gam(data=data[kfold!=fold & kfold!="Test" & BadData==F,],
                            formula = formula,
                            weights = w,...)
        if(model_res2){
          gams[[paste0(fold,"_r")]] <- gam(data=cbind(data[kfold!=fold & kfold!="Test" & BadData==F,],data.table(gam_res2=residuals(gams[[fold]])^2)),
                                           formula = formula_res2,
                                           weights = w,...)
        }
      }
      
      ## Predictions
      data$gam_pred[data$kfold==fold] <- predict(gams[[fold]],newdata = data[kfold==fold,])  
      
    }
  }else{
    warning("Skipping gam modelling and using provided gam_pred.")
    gams <- list()
  }
  
  
  data[,gam_resid:=get(as.character(formula[[2]]))-gam_pred]
  
  ## Container for predictive quantiles
  predqs <- data.frame(matrix(NA,ncol = length(quantiles), nrow = nrow(data)))
  colnames(predqs) <- paste0("q",100*quantiles)
  
  
  ## Qunatile regression
  if(!is.null(formula_qr)){
    ## Use user-specified model equation for quantile regression
    formula_qr <- reformulate(attr(terms(formula_qr),"term.labels"),response = "gam_resid")
    for(fold in unique(data$kfold)){
      print(paste0("MQR, kfold=",fold))
      for(i in 1:length(quantiles)){
        
        lqr_fit <- rq(formula_qr,
                      tau = quantiles[i],
                      data = data[kfold!=fold & kfold!="Test" & BadData==F,],
                      weights = data[kfold!=fold & kfold!="Test" & BadData==F,w],
                      method = "br")
        
        predqs[data$kfold==fold,i] <- data[kfold==fold,gam_pred] + predict.rq(lqr_fit,data[kfold==fold,])
      }
    }
  }else{
    ## QR with features from GAM
    for(fold in unique(data$kfold)){
      print(paste0("MQR, kfold=",fold))
      ## Get training Data
      train <- predict(gams[[fold]],newdata = data[kfold!=fold & kfold!="Test" & BadData==F,],type = "terms")
      if(model_res2){
        train2 <- predict(gams[[paste0(fold,"_r")]],newdata = data[kfold!=fold & kfold!="Test" & BadData==F,],type = "terms")
        # Only need to retrun smooth terms as linear terms included already...
        train2 <- train2[,grep("\\(",colnames(train2)),drop=F]
        colnames(train2) <- paste0(colnames(train2),"_r")
        train <- cbind(train,train2); rm(train2)
      }
      train <- cbind(as.data.table(train),data[kfold!=fold & kfold!="Test" & BadData==F,.(gam_resid)])
      
      ## Out-of-sample data
      test_cv <- as.data.table(predict(gams[[fold]],newdata = data[kfold==fold,],type = "terms"))
      if(model_res2){
        test_cv2 <- predict(gams[[paste0(fold,"_r")]],newdata = data[kfold==fold,],type = "terms")
        test_cv2 <- test_cv2[,grep("\\(",colnames(test_cv2)),drop=F]
        colnames(test_cv2) <- paste0(colnames(test_cv2),"_r")
        test_cv <- cbind(test_cv,test_cv2); rm(test_cv2)
      }
      
      for(i in 1:length(quantiles)){
        ## Fit QR model
        lqr_fit <- rq(formula = gam_resid ~ .,
                      tau = quantiles[i],
                      data = train,
                      weights = data[kfold!=fold & kfold!="Test" & BadData==F,w],
                      method = "br")
        
        
        
        ## Make predictions
        predqs[data$kfold==fold,i] <- data[kfold==fold,gam_pred] + predict.rq(lqr_fit,test_cv)
        
      }
    }
  }
  
  class(predqs) <- c("MultiQR","data.frame")
  
  
  if(Sort){
    predqs <- SortQuantiles(data = predqs,Limits = SortLimits)
  }
  
  return(list(predqs=predqs,
              gams=gams,
              gam_pred=data$gam_pred))
  
}






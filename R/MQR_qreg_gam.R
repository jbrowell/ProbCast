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
#' numbered/labeled folds and "Test" for test data.
#' @param formala A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right passed to \code{gam()}
#' or \code{bam()} from \code{mgcv}.
#' @param formula_qr Formula for linear quantile regression model for GAM residuals. Term \code{gam_pred}
#' is the prediction from the above GAM may be included in this formula.
#' @param model_res2 If \code{TRUE} also model squared residuals of GAM using a GAM. Defaults to \code{FALSE}.
#' @param formula_res2 Formula for GAM to predict squared residuals.
#' @param quantiles The quantiles to fit models for.
#' @param gbm_params List of parameters to be passed to \code{fit.gbm()}.
#' @param CVfolds Control for cross-validation if not supplied in \code{data}.
#' @param Sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param SortLimits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @details Returns a \code{list} comprising predictive quantiles, GAM models, and deterministic predictions
#' from GAMs.
#' 
#' The returned predictive quantiles and GAM predictions are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' @return Quantile forecasts in a \code{MultiQR} object.
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

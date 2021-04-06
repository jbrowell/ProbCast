###############
## GAM+MQR ####
###############
#' Multiple Quantile Regression Using Generalised Additive Models and Linear Quantile Regression
#'
#' This function fits multiple conditional linear quantile regression models to the residuals of
#' a generalised additive model using \code{mgcv} and with facilities for cross-validation.
#' Quantile regression may be performed using user-specified formula or the design matrix of
#' the fitted GAM.
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
#' @param cv_folds Control for cross-validation with various options, either:
#' \itemize{
#'  \item the column name of the fold index supplied in data. Observations and inputs 
#'   in the index labeled "Test" will serve as test data and held out in model training.
#'  \item an integer giving the number of cross validation folds to generate. Folds are constructed as block chunks. 
#'  Default behaviour is 5 folds.
#'  \item NULL indicates that no cross validation should be performed and the returned model is trained on all \code{data}.
#' }
#' @param use_bam If \code{TRUE} (default) then GAM is fit using (\code{bam()}) in stead
#' of \code{gam()}. \code{bam} is better suited to large datasets but not all
#' \code{gam} model options are available with \code{bam}. Alternative smooths, such as
#' cubic regression splines aid faster estimation than the default smooth \code{bs="tp"}.
#' See \code{bam()} documentation for further details.
#' @param exclude_train A column name in \code{data} indicating if a row should be excluded from model
#' training, i.e. if it contains bad data (will be coerced to \code{logical}). Alterntively,
#' an \code{integer} or \cpde{logical} vector with length equal to the number of rows in \code{data} indicating
#' the same. Rows labeled \code{TRUE} are excluded from model training.
#' @param sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param ... Additional agruments passter to \code{gam()} (or \code{bam()}).
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
                     exclude_train=NULL,
                     # w=rep(1,nrow(data)), # Not working...
                     sort=T,sort_limits=NULL,
                     ...){
  
  ## To do:
  # - Add use of "cluster" option
  # - Add warning if bam() used with default s(bs="tp"), this slows things alot.
  
  # Check compatible options:
  if(!is.null(formula_qr) & model_res2){
    stop("Additinal GAM effects for modelling squared residuals only available if quantile regression is based on effects of main GAM, i.e. \"formula_qr=NULL\"")}
  if(!model_res2 & !is.null(formula_res2)){
    stop("formula_res2 not required as model_res2==F")
  }
  # set-up cv folds & do checks
  cv_labs <- cv_control(data = data,cv_folds = cv_folds)
  
  # exclude points from training? & do checks
  exclude_idx <- exclude_fun(data = data,exclude_train = exclude_train)
  
  
  # Special names...
  if("gam_resid" %in% names(data)){warning("data$gam_resid will be overwritten!")}
  if("gam_pred" %in% names(data)){warning("data$gam_pred will be overwritten!")}
  
  ## GAM for conditional expectation (and possible squared residuals)
  
  FINAL_OUTPUT <- list(mqr_pred=NULL,
                       call=match.call(),
                       kfold_index = cv_labs$idx,
                       model_names = cv_labs$fold_loop,
                       exclude_index = exclude_idx,
                       models=list(gam_pred=data.table::copy(data[,.(y=get(as.character(formula[[2]])))]),
                                   gams = list(),
                                   rqs = list()),
                       sorted=list(sort=sort,
                                   sort_limits=sort_limits),
                       default_model=if("Test"%in%cv_labs$fold_loop){"Test"}else{cv_labs$fold_loop[1]})
  
  class(FINAL_OUTPUT) <- c("qreg_gam",class(FINAL_OUTPUT))
  
  formula_res2 <- reformulate(attr(terms(formula_res2),"term.labels"),response = "gam_res2")
  
  for(fold in FINAL_OUTPUT$model_names){
    
    print(paste0("GAM, kfold=",fold))
    
    ## Fit gams using either gam() or bam()
    gam_fit_method <- if(use_bam){bam}else{gam}
    
    # GAM for conditional expectation
    FINAL_OUTPUT$models$gams[[fold]] <- gam_fit_method(data=data[FINAL_OUTPUT$kfold_index!=fold & FINAL_OUTPUT$kfold_index!="Test" & FINAL_OUTPUT$exclude_index==0,],
                                                       formula = formula,...)
    # GAM for squared residuals
    if(model_res2){ 
      temp_gam_res2 <- data.table(gam_res2=(FINAL_OUTPUT$models$gam_pred[FINAL_OUTPUT$kfold_index!=fold & FINAL_OUTPUT$kfold_index!="Test" & FINAL_OUTPUT$exclude_index==0,y] - 
                                              predict(FINAL_OUTPUT$models$gams[[fold]],
                                                      newdata = data[FINAL_OUTPUT$kfold_index!=fold & FINAL_OUTPUT$kfold_index!="Test" & FINAL_OUTPUT$exclude_index==0,]))^2)  
      
      FINAL_OUTPUT$models$gams[[paste0(fold,"_r")]] <- gam_fit_method(data=cbind(data[FINAL_OUTPUT$kfold_index!=fold & FINAL_OUTPUT$kfold_index!="Test" & FINAL_OUTPUT$exclude_index==0,],temp_gam_res2),
                                                                      formula = formula_res2,...)
      rm(temp_gam_res2)
    }
    
    
    ## Out-of-sample cross-validation predictions
    if(is.null(cv_folds)){
      FINAL_OUTPUT$models$gam_pred[,gam_pred:=predict(FINAL_OUTPUT$models$gams[[fold]],newdata = data)]  
    }else{
      FINAL_OUTPUT$models$gam_pred[FINAL_OUTPUT$kfold_index==fold,
                                   gam_pred:=predict(FINAL_OUTPUT$models$gams[[fold]],newdata = data[FINAL_OUTPUT$kfold_index==fold,])]  
    }
  }
  
  ## Out-of-sample cross-validation residuals
  FINAL_OUTPUT$models$gam_pred[,gam_resid:=y-gam_pred]
  
  ## Quantile regression ##
  FINAL_OUTPUT <- qreg_gam.add_quantiles(object=FINAL_OUTPUT,
                                         data=data,
                                         quantiles = quantiles)
  
  return(FINAL_OUTPUT)
  
}



#' Add new qunatile regression models and predictions to \code{qreg_gam} model.
#'
#' This function adds new conditional quantile models and associated predictions
#' to an exisiting \code{qreg_gam} object. This is a much faster approach to add
#' new quantiles than re-estimating the main GAM using \code{qreg_gam()}. This
#' function is also called within \code{qreg_gam()}.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param object An \code{qreg_gam} object.
#' @param data The data used to fit \code{object}
#' @param quantiles The new quantiles to be added to \code{object}. Models for quantiles already
#' in \code{object} will not be re-estimated and a warning will be thrown.
#' @details Adds new qunatile models to \code{qreg_gam} object. 
#' @return An updated \code{qreg_gam}.
#' @export
qreg_gam.add_quantiles <- function(object, data, quantiles){
  
  
  ## Check quantiles don't already exisit
  if(any(quantiles %in% (as.numeric(gsub("q","",colnames(object$mqr_pred)))/100))){
    
    quantiles <- quantiles[-which(quantiles %in% (as.numeric(gsub("q","",colnames(object$mqr_pred)))/100))]
    
    if(length(quantiles)==0){
      warning("All quantiles already exisit are are not re-estimated.")
      return(object)
    }
    else{
      warning("Some quantiles already exisit and are not re-estimated.")
    }
  }
  
  
  ## Need to assign gam_resid for modelling
  if("gam_resid" %in% names(data)){warning("data$gam_resid will be overwritten!")}
  if("gam_pred" %in% names(data)){warning("data$gam_pred will be overwritten!")}
  
  data$gam_resid <- object$models$gam_pred$gam_resid
  data$gam_pred <- object$models$gam_pred$gam_pred
  
  
  ## Container for predictive quantiles
  predqs <- data.table(matrix(as.numeric(NA),ncol = length(quantiles), nrow = nrow(data)))
  colnames(predqs) <- paste0("q",100*quantiles)
  
  if(!is.null(object$call$formula_qr)){
    ## Use user-specified model equation for quantile regression
    formula_qr <- eval(object$call$formula_qr,envir = environment())
    formula_qr <- reformulate(attr(terms(as.formula(formula_qr)),"term.labels"),response = "gam_resid")
    
    for(fold in object$model_names){
      print(paste0("MQR, kfold=",fold))
      for(i in 1:length(quantiles)){
        
        rq_model <- rq(formula_qr,
                       tau = quantiles[i],
                       data = data[object$kfold_index!=fold & object$kfold_index!="Test" & object$exclude_index==0,],
                       method = "br")
        
        ## Throw away unneeded entries in "rq" as take up a lot of space!!!
        rq_model <- rq_model[c("coefficients","terms","xlevels","contrasts")]
        class(rq_model) <- "rq"
        object$models$rqs[[fold]][[paste0("q",100*quantiles[i])]] <- rq_model
        
        if(is.null(object$call$cv_folds)){
          predqs[,i] <- object$models$gam_pred[,gam_pred] +
            predict.rq(object$models$rqs[[fold]][[paste0("q",100*quantiles[i])]],data)
        }else{
          predqs[object$kfold_index==fold,i] <- object$models$gam_pred[object$kfold_index==fold,gam_pred] +
            predict.rq(object$models$rqs[[fold]][[paste0("q",100*quantiles[i])]],data[object$kfold_index==fold,])
        }
      }
    }
  }else{
    ## QR with features from GAM
    for(fold in object$model_names){
      print(paste0("MQR, kfold=",fold))
      ## Get training Data
      train <- predict(object$models$gams[[fold]],newdata = data[object$kfold_index!=fold & object$kfold_index!="Test" & object$exclude_index==0,],type = "terms")
      if(!is.null(object$call$formula_res2)){
        train2 <- predict(object$models$gams[[paste0(fold,"_r")]],newdata = data[object$kfold_index!=fold & object$kfold_index!="Test" & object$exclude_index==0,],type = "terms")
        # Only need to retrun smooth terms as linear terms included already...
        train2 <- train2[,grep("\\(",colnames(train2)),drop=F]
        if(ncol(train2)>0){
          colnames(train2) <- paste0(colnames(train2),"_r")
          train <- cbind(train,train2)
        }
        rm(train2)
      }
      train <- cbind(data.table(train),data[object$kfold_index!=fold & object$kfold_index!="Test" & object$exclude_index==0,.(gam_resid)])
      
      ## Out-of-sample data
      
      test_cv <- predict(object$models$gams[[fold]],
                         newdata = if(is.null(object$call$cv_folds)){data}else{data[object$kfold_index==fold,]},
                         type = "terms")
      if(!is.null(object$call$formula_res2)){
        test_cv2 <- predict(object$models$gams[[paste0(fold,"_r")]],
                            newdata = if(is.null(object$call$cv_folds)){data}else{data[object$kfold_index==fold,]},
                            type = "terms")
        test_cv2 <- test_cv2[,grep("\\(",colnames(test_cv2)),drop=F]
        if(ncol(test_cv2)>0){
          colnames(test_cv2) <- paste0(colnames(test_cv2),"_r")
          test_cv <- cbind(test_cv,test_cv2)
        }
        rm(test_cv2)
      }
      test_cv <- data.table(test_cv)
      
      
      for(i in 1:length(quantiles)){
        ## Fit QR model
        object$models$rqs[[fold]][[paste0("q",100*quantiles[i])]] <- 
          rq(formula = gam_resid ~ .,
             tau = quantiles[i],
             data = train,
             method = "br")
        
        ## Make predictions
        if(is.null(object$call$cv_folds)){
          predqs[,i] <- object$models$gam_pred[,gam_pred] +
            predict.rq(object$models$rqs[[fold]][[paste0("q",100*quantiles[i])]],
                       newdata=test_cv)
        }else{
          predqs[object$kfold_index==fold,i] <- object$models$gam_pred[object$kfold_index==fold,gam_pred] +
            predict.rq(object$models$rqs[[fold]][[paste0("q",100*quantiles[i])]],
                       newdata=test_cv)
        }
      }
    }
  }
  
  
  ## Combine with existing quantiles
  object$mqr_pred <- cbind(object$mqr_pred,predqs)
  
  ## Order columns
  object$mqr_pred <- object$mqr_pred[,order(as.numeric(gsub("q","",colnames(object$mqr_pred)))),with=F]
  
  ## Assigne MultiQR class
  if(class(object$mqr_pred)[1]!="MultiQR"){
    class(object$mqr_pred) <- c("MultiQR",class(object$mqr_pred))
  }
  
  ## Sort quantiles
  if(object$sorted$sort){
    object$mqr_pred <- SortQuantiles(data = object$mqr_pred,
                                     Limits = object$sorted$sort_limits)
  }
  
  data$gam_resid <- NULL
  data$gam_pred <- NULL
  
  
  return(object)
}


#' Predict from model based on Generalised Additive Model and Linear Quantile Regression
#'
#' This function predicts from multiple conditional linear quantile regression models of the residuals of
#' a generalised additive model fit using \code{qreg_gam}.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param object An \code{qreg_gam} object containing the model to predict from.
#' @param newdata A data frame or data table containing the values of the model
#' covariates at which predictions are required. 
#' @param quantiles The probability levels at which quantile predictions should
#' be produced.
#' @param model_name The name of the model in \code{object} to be used for prediction. E.g.
#' Specific cross-vlaidation fold, test model or some other version.
#' @param sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. 
#' Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @details Predict method for multiple quantile regression models of the class \code{qreg_gam}. 
#' @return A list with elements \code{gam_pred}, deterministic predictions (conditional expectation) from
#' main GAM model, and \code{mqr_pred}, multiple predictive quantiles in a \code{MultiQR} object.
#' @export
predict.qreg_gam <- function(object,
                             newdata = NULL,
                             quantiles = NULL,
                             model_name=NULL,
                             sort = T,
                             sort_limits = NULL){
  
  ## Checks
  # Class
  if(class(object)[1]!="qreg_gam"){stop("object of wrong class, expecting \"qreg_gam\"")}
  # Correct columns in newdata
  # Availability of quantile models
  
  ## Use default modeul unless specified
  if(is.null(model_name)){
    model_name <- object$default_model
  }
  
  if(is.null(quantiles)){
    quantiles <- sort(as.numeric(gsub(pattern = "q",replacement = "",names(object$models$rqs[[model_name]])))/100)
  }
  if(!all(quantiles %in% (as.numeric(gsub(pattern = "q",replacement = "",names(object$models$rqs[[model_name]])))/100))){
    stop("Models not availalbe for all requested quantiles.")
  }
  
  ## Initialise containers
  predqs <- data.table(matrix(as.numeric(NA),ncol = length(quantiles), nrow = nrow(newdata)))
  colnames(predqs) <- paste0("q",100*quantiles)
  
  OUTPUT <- list(gam_pred=predict(object$models$gams[[model_name]],newdata = newdata),
                 mqr_pred=NULL)
  newdata[,gam_pred:=OUTPUT$gam_pred]
  
  ## Qunatile regression
  if(!is.null(object$call$formula_qr)){
    ## Use user-specified model equation for quantile regression
    for(i in 1:length(quantiles)){
      predqs[,i] <- OUTPUT$gam_pred +
        predict.rq(object = object$models$rqs[[model_name]][[paste0("q",100*quantiles[i])]],
                   newdata = newdata)
    }
  }else{
    ## QR with features from GAM
    newdata_terms <- predict(object$models$gams[[model_name]],newdata = newdata,type = "terms")
    if(!is.null(object$call$formula_res2)){
      newdata_terms2 <- predict(object$models$gams[[paste0(model_name,"_r")]],newdata = newdata,type = "terms")
      newdata_terms2 <- newdata_terms2[,grep("\\(",colnames(newdata_terms2)),drop=F]
      colnames(newdata_terms2) <- paste0(colnames(newdata_terms2),"_r")
      newdata_terms <- cbind(newdata_terms,newdata_terms2); rm(newdata_terms2)
    }
    newdata_terms <- data.table(newdata_terms)
    
    ## Make predictions
    for(i in 1:length(quantiles)){
      predqs[,i] <- OUTPUT$gam_pred +
        predict.rq(object = object$models$rqs[[model_name]][[paste0("q",100*quantiles[i])]],
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





#' Update a GAM with new training data.
#'
#' This function updates a \code{qreg_gam} model with new training data.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param object An \code{qreg_gam} object containing the model to update.
#' @param newdata A data frame or data table containing the values of the model
#' covariates and target variable required to update the model. 
#' @param model_name The name of the model in \code{object} to be updated.
#' @details Update the \code{bam} component of an \code{qreg_gam} model. 
#' @return An updated \code{qreg_gam} model.
#' @export
qreg_gam.update <- function(object,newdata,model_name=NULL){
  
  
  ## Check class of object
  if(class(object)[1]!="qreg_gam"){stop("object of wrong class, expecting \"qreg_gam\"")}
  
  ## Use default model unless specified
  if(is.null(model_name)){
    model_name <- object$default_model
  }
  
  ## Check GAMs were fit using bam()
  if(!"bam" %in% class(object$models$gams[[model_name]])){
    stop("Only suitable for bam(). Select \"use_bam=T\" when fitting.")
  }
  
  
  
  ## Check model can be updated (from bam.update with extra explanation)
  if(is.null(object$models$gams[[model_name]]$qrx)){
    stop("Model can not be updated. Possible reason: bam option descrete=T, must be =F for updating.")
  }
  
  object$models$gams[[model_name]] <- bam.update(b = object$models$gams[[model_name]],
                                                 data = newdata)
  
  return(object)
  
} 


#' Print \code{qreg_gbm} object
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param x object of class \code{qreg_gam} obtained from the function \code{qreg_gbm()}.
#' @details this function prints details of the call to \code{qreg_gam()},
#' @keywords Quantile Regression
#' @method print qreg_gam
#' @export    
print.qreg_gam <- function(x, ...){
  cat("\n")
  cat("\t Multiple Quantile regression via GBM \n")
  cat("\n")
  cat("Call:\n", deparse(x$call), "\n\n", sep = "")
  
  cat("\n")
  cat("Names of models fitted: ", paste0(x$model_names,collapse = ", "),"\n")
  cat("Default prediction model: ", x$default_model,"\n")
  cat("Quantiles fitted for each model: ", paste0(names(x$mqr_pred),collapse = ", "),"\n")
  cat("\n")
  invisible(x)
}


#' Summary \code{qreg_gbm} object
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param x object of class \code{qreg_gam} obtained from the function \code{qreg_gbm()}.
#' @details this function gives a more detailed summary of, \code{x}, a \code{qreg_gam} object,
#' than \code{print()}
#' @keywords Quantile Regression
#' @method summary qreg_gam
#' @export
summary.qreg_gam <- function(x, ...){
  
  print(x)
  
  cat("\n")
  
  if(!is.null(x$mqr_pred)){
    cat("Out-of-sample kfold cross validation predictions:\n")
    print(data.table::data.table(x$mqr_pred))
    
    cat("\n")
    cat("class of kfold predictions: ", class(x$mqr_pred), "\n")
    cat("sorted quantile predictions?: ", x$sorted$sort, "\n")
    if(x$sorted$sort){cat("sort limits: ", deparse(x$sorted$sort_limits), "\n")}
    
  }
  
  cat("\n-----\n\n")
  cat("Summary of default gam model:")
  cat("\n")
  print(summary(x$models$gams[[x$default_model]]))
  cat("\n")
  cat("For more diagnostics see gam.check()")
  cat("\n-----\n")
  
}

## DEPRECIATED ####
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
  
  warning("This function is depreciated and has been replaced by qreg_gam().")
  
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






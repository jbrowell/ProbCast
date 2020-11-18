#' Multiple Quantile Regression Using Gradient Boosted Decision Trees
#'
#' This function fits multiple boosted quantile regression trees 
#' using \code{gbm} with facilities for cross-validation.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}; Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory
#' variables. May optionally contain a column called "kfold" with
#' numbered/labeled folds and "Test" for test data.
#' @param formala A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right
#' @param quantiles The quantiles to fit models for.
#' @param cv_folds Control for cross-validation if not supplied in \code{data}.
#' @param perf.plot Plot GBM performance?
#' @param pred_ntree predict using a user-specified tree.
#' If unspecified an out-of-the bag estimate will be used unless internal
#' gbm cross-validation folds are specified in \code{...}.
#' @param cores the number of available cores. Defaults to one, i.e. no parallelisation, although in this case the user
#'  must still specify \code{pckgs} if applicable.
#' @param pckgs specify additional packages required for
#' each worker (e.g. c("data.table") if data stored as such).
#' @param sort Sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and 
#' lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param save_models_path Path to save models. Model details and file extension pasted onto this string. 
#' Defaults to \code{NULL}, i.e. no model save.
#' @param return_op_model return an operational model trained on \strong{all available data} for prediction? see \code{predict.pc_gbm()}.
#' @param ... Additional arguments passed to \code{gbm()}.
#' @details The returned predictive quantiles are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data. 
#' 
#' When \code{return_op_model = TRUE} the function returns a list containing the \code{MultiQR} object
#' and an operational model for generating predictions.
#' @return Quantile forecasts in a \code{MultiQR} object. Optionally a list containing a \code{MultiQR}
#' object and a model for operational prediction
#' @keywords Quantile Regression
#' @importFrom foreach %dopar%
#' @export
mqr_qreg_gbm <- function(data,
                         formula,
                         quantiles = c(0.25,0.5,0.75),
                         cv_folds = NULL,
                         perf.plot = F,
                         pred_ntree = NULL,
                         cores = 1,
                         pckgs = NULL,
                         sort = T,
                         sort_limits = NULL,
                         save_models_path = NULL,
                         return_op_model = FALSE,
                         ...){
  
  # to-do
  ## clear up auto-cv --- cv_folds = NULL and no-kfold in data? convert to function
  ### no kfold capability?
  #### issue/target/fold indexed mqr object?
  
  # Set-up Cross-validation
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
  
  
  # model folds check -> as.character for op model option
  if(is.numeric(data$kfold)){
    data$kfold <- as.character(data$kfold)
  }
  
  # model training folds
  if(return_op_model){
    fold_mods <- c(unique(data$kfold),"op_model")
  } else{
    fold_mods <- unique(data$kfold)
  }
  
  
  # set up parallel workers, defaults to one worker....
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  #set up progress bar
  iterations <- length(quantiles)
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  gc()
  
  # fit the models returns list quantile --> kfold
  q_mods <- foreach::foreach(q = quantiles,.packages = c("gbm",pckgs),.options.snow = opts) %dopar% {
    
    temp_gbm <- list()
    
    for(fold in fold_mods){
      
      if(fold=="op_model"){
        
        # need to define as object for save option
        temp <- gbm::gbm(formula=formula,
                         data = data[!is.na(data[[formula[[2]]]]),],
                         distribution = list(name="quantile",alpha=q),
                         ...)
        
        
      } else{
        
        temp <- gbm::gbm(formula=formula,
                         data = data[data$kfold!=fold & data$kfold!="Test" & !is.na(data[[formula[[2]]]]),],
                         distribution = list(name="quantile",alpha=q),
                         ...)
      }
      
      
      # save model?
      if(!is.null(save_models_path)){
        if(fold=="op_model"){
          try(save(temp,file = paste0(save_models_path,"_q",100*q,"_",fold,".rda")))}
        else{
          try(save(temp,file = paste0(save_models_path,"_q",100*q,"_fold",fold,".rda")))
        }
      }
      
      temp_gbm[[fold]] <- temp

      
      
    }
    
    return(temp_gbm)
    
  }
  
  close(pb)
  parallel::stopCluster(cl)
  names(q_mods) <- paste0("q",100*quantiles)
  
  
  # re-arrange structure so  kfold --> quantiles
  # assign new class to each kfold sub model
  q_mods <- lapply(fold_mods,function(x){
    foldmod <- lapply(q_mods,function(z){
      z[[x]]
    })
    class(foldmod) <- "pc_gbm"
    return(foldmod)
  })
  names(q_mods) <- fold_mods
  


  ### create container for cv output
  predqs <- data.frame(matrix(NA,ncol = length(quantiles), nrow = nrow(data)))
  colnames(predqs) <- paste0("q",100*quantiles)
  
  for(fold in unique(data$kfold)){
    
    predqs[data$kfold==fold,] <- predict(q_mods[[fold]],
                                         newdata = data[data$kfold==fold,],
                                         quantiles = NULL,
                                         pred_ntree = pred_ntree,
                                         perf.plot = perf.plot,
                                         sort = sort,
                                         sort_limits = sort_limits)
    
  }
  
  # class of new cv object
  class(predqs) <- c("MultiQR","data.frame")
 
  
  if(return_op_model){
    return(list(cv_mqr = predqs,
                op_mod = q_mods[["op_model"]]))
  } else{
    return(predqs)
  }
  
  
}



#' Multiple Quantile Regression Using Gradient Boosted Decision Trees (depreciated)
#'
#' This function is now depreciated and may be removed in future versions of this package.
#' Use \code{mqr_qreg_gbm()} instead.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}; Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory
#' variables. May optionally contain a collumn called "kfold" with
#' numbered/labeled folds and "Test" for test data.
#' @param formala A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right
#' @param quantiles The quantiles to fit models for.
#' @param gbm_params List of parameters to be passed to \code{fit.gbm()}.
#' @param CVfolds Control for cross-validation if not supplied in \code{data}.
#' @param perf.plot Plot GBM performance?
#' @param pred_ntree predict using a user-specified tree.
#' If unspecified an out-of-the bag estimate will be used unless interval
#' gbm cross-validation folds are specified in \code{gbm_params}.
#' @param parallel \code{boolean} parallelize cross-validation process?
#' Parallelisation is over cross-validation folds by default, optionally over
#' quantiles by setting code{para_over_q=T}.
#' @param pckgs if \code{parallel} is TRUE then  specify packages required for
#' each worker (e.g. c("data.table) if data stored as such).
#' @param cores if \code{parallel} is TRUE then number of available cores
#' @param para_over_q if \code{parallel} is TRUE then paralellize over quantiles?
#' Defalts to FALSE i.e."kfold".
#' @param Sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param SortLimits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @details The returned predictive quantiles are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' @return Quantile forecasts in a \code{MultiQR} object.
#' @keywords Quantile Regression
#' @export
MQR_gbm <- function(data,
                    formula,
                    quantiles=c(0.25,0.5,0.75),
                    CVfolds=NULL,
                    gbm_params=list(...),
                    perf.plot=F,
                    parallel = F,
                    pred_ntree = NULL,
                    cores = NULL,
                    pckgs = NULL,
                    para_over_q = FALSE,
                    Sort=T,
                    SortLimits = NULL){
  
  
  warning("function depreciated and may be removed in future updates. Update to ProbCast::mqr_qreg_gbm()")
  
  
  if(parallel){
    
    if(isFALSE(para_over_q)){
      
      warning("parallelism is now only spported over quantiles. Calculating over quantiles...")
      
    }
  
  
    do.call(mqr_qreg_gbm,c(list(data = data, 
                                formula = formula,
                                quantiles=quantiles,
                                cv_folds=CVfolds,
                                perf.plot=perf.plot,
                                pred_ntree = pred_ntree,
                                cores = cores,
                                pckgs = pckgs,
                                sort = Sort,
                                sort_limits= SortLimits,
                                save_models_path=NULL),
                           gbm_params))
    
    
  } else{
    
    
    do.call(mqr_qreg_gbm,c(list(data = data, 
                                formula = formula,
                                quantiles=quantiles,
                                cv_folds=CVfolds,
                                perf.plot=perf.plot,
                                pred_ntree = pred_ntree,
                                cores = 1,
                                pckgs = pckgs,
                                sort = Sort,
                                sort_limits= SortLimits,
                                save_models_path=NULL),
                           gbm_params))
    
    
    
  }
    
  
  
}




#' Predict method for Multiple Quantile Regression from Gradient Boosted Decision Trees.
#'
#' This function returns multiple quantile predictions as an object of class \code{MultiQR}
#' based on a ProbCast gbm fitted model. S3 Method for for \code{pc_gbm}
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param object 	object of class \code{pc_gbm} obtained from the function \code{mqr_qreg_gbm()}
#' @param newdata data.frame of observations for which to make predictions
#' @param quantiles The quantiles to predict. Default is all the quantiles present in \code{object}
#' @param pred_ntree predict using a user-specified tree.
#' If unspecified an out-of-the bag estimate will be used unless internal
#' gbm cross-validation folds are specified in \code{mqr_qreg_gbm()}
#' @param perf.plot plot GBM performance if \code{pred_ntree = NULL}?
#' @param sort sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and 
#' lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param ... additional arguments; not currently used.
#' @details this function returns predictive quantiles for each row in \code{newdata},
#' the result is returned as a \code{MultiQR} object 
#' @return Quantile forecasts in a \code{MultiQR} object.
#' @keywords Quantile Regression
#' @export
predict.pc_gbm <- function(object,
                           newdata = NULL,
                           quantiles = NULL,
                           pred_ntree = NULL,
                           perf.plot = FALSE,
                           sort = T,
                           sort_limits = NULL,
                           ...) {
  
  
  if(is.null(quantiles)){
    quantiles <- names(object)
  }
  
  if(is.numeric(quantiles)){
    quantiles <- paste0("q",quantiles*100)
  }
  
  
  pred <- data.frame(sapply(quantiles,function(q){
    
    ### Save out-of-sample predictions
    if(is.null(pred_ntree)){
      gbm::predict.gbm(object[[q]],
                       newdata = newdata,
                       n.trees = gbm::gbm.perf(object[[q]],plot.it = perf.plot))
    } else{
      
      gbm::predict.gbm(object[[q]],
                       newdata = newdata,
                       n.trees = pred_ntree)
      
    }
    
    
  }))
  
  class(pred) <- c("MultiQR","data.frame")
  
  
  if((sort) & (length(quantiles) != 1)){
    pred <- SortQuantiles(data = pred,Limits = sort_limits)
  }
  
  
  return(pred)
  
  
  
}







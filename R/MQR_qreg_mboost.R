#' Multiple Quantile Regression using \code{mboost}
#'
#' This function fits multiple quantile regreesion models using \code{mboost}, with
#' facilities for cross-validation. \code{mboost} accommodates both generalised
#' additive models, decision trees and other learners. See \code{?mboost} for more
#' details.
#' 
#' Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' 
#' @param data A \code{data.frame} containing target and explanatory variables.
#' May optionally contain a column called "kfold" with numbered/labeled folds
#' and "Test" for test data.
#' @param formaula A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right.
#' @param quantiles A vector with length>=2 containing the quantiles to fit models for.
#' @param cv_folds Control for cross-validation if not supplied in \code{data}.
#' @param bc_mstop for boosting control, the number of boosting iterations
#' @param bc_nu for boosting control, the shrinkage or step size with a value between 0 and 1
#' @param w an optional numeric vector of weights to be used in the fitting process.
#' @param cores the number of available cores. Defaults to one, i.e. no parallelisation, although in 
#' this case the user must still specify \code{pckgs} if applicable.
#' @param pckgs specify additional packages required for
#' each worker (e.g. c("data.table") if data stored as such).
#' @param sort Sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles 
#' to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param save_models_path Path to save models. Model details and file extension pasted onto this string.
#' Defaults to \code{NULL}, i.e. no model save.
#' @param ... extra hyper-parameters to be passed to \code{mboost()}.
#' @details The returned predictive quantiles are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' @return Quantile forecasts in a \code{MultiQR} object.
#' @keywords Quantile Regression
#' @importFrom foreach %dopar%
#' @export
mqr_qreg_mboost <- function(data,
                            formula,
                            quantiles=c(0.25,0.5,0.75),
                            cv_folds=NULL,
                            bc_mstop=100,
                            bc_nu=0.1,
                            w=rep(1,nrow(data)),
                            cores = 1,
                            pckgs = NULL,
                            sort=T,
                            sort_limits=NULL,
                            save_models_path=NULL,
                            ...){
  
  ## to-do
  ## -- clear up auto-cv
  ## -- rolling window
  ## -- data.table, dependency with BadData?
  ## -- boosting control?
  ## -- check single quantiles warning
  
  if(length(quantiles)<2){
    stop("length(quantiles)<2. This function takes issue with single quantiles.")
  }
  
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
  
  # discuss with jethro for other functions as well...
  if(!"BadData" %in% colnames(data)){
    ### assumed data.table?
    # data[,BadData:=F]
    data$BadData <- rep(FALSE,nrow(data))
  }
  
  
  ### create container for output
  predqs <- data.frame(matrix(NA,ncol = length(quantiles), nrow = nrow(data)))
  colnames(predqs) <- paste0("q",100*quantiles)
  
  
  ### set up parallel workers, defaults to one worker....
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  #set up progress bar
  iterations <- length(quantiles)
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  gc()
  
  
  
  ### fit the models
  qpred <- foreach::foreach(q = quantiles,.packages = c("mboost",pckgs),.options.snow = opts) %dopar% {
    
    
    
    pred <- list()
    
    for(fold in unique(data$kfold)){
      
      
      temp_model <- mboost::mboost(formula=formula,
                                   data=data[data$kfold!=fold & data$kfold!="Test" & data$BadData==F & !is.na(data[[formula[[2]]]]),],
                                   family = mboost::QuantReg(tau=q),
                                   control= mboost::boost_control(mstop = bc_mstop,
                                                                  nu = bc_nu),
                                   weights=w[data$kfold!=fold & data$kfold!="Test" & data$BadData==F & !is.na(data[[formula[[2]]]])],
                                   ...)
      
      if(!is.null(save_models_path)){
        try(save(temp_model,file = paste0(save_models_path,"_q",100*q,"_fold",fold,".rda")))
      }
      
      
      
      pred[[fold]] <- mboost::predict.mboost(temp_model,
                                             newdata = data[data$kfold==fold,])
      
      
    }
    
    
    return(pred)
    
    
  }
  
  
  close(pb)
  parallel::stopCluster(cl)
  
  names(qpred) <- paste0("q",100*quantiles)
  
  for(q in quantiles){
    for(fold in unique(data$kfold)){# Loop over CV folds and test data
      predqs[[paste0("q",100*q)]][data$kfold==fold] <- qpred[[paste0("q",100*q)]][[fold]]
      
    }
  }
  
  
  
  class(predqs) <- c("MultiQR","data.frame")
  
  
  if((sort) & (length(quantiles) != 1)){# if only one quantile specified, sortquantiles transposes predqs (also no need to sort).
    predqs <- SortQuantiles(data = predqs,Limits = sort_limits)
  }
  
  
  
  return(predqs)
  
}























#' Multiple Quantile Regression using \code{mboost} (depreciated)
#'
#' This function is now depreciated and may be removed in future versions of this package.
#' Use \code{mqr_qreg_mboost()} instead.
#' 
#' Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' 
#' @param data A \code{data.frame} containing target and explanatory variables.
#' May optionally contain a column called "kfold" with numbered/labeled folds
#' and "Test" for test data.
#' @param formaula A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right.
#' @param quantiles A vector with length>=2 containing the quantiles to fit models for.
#' @param model_params List of parameters to be passed to \code{fit.gbm()}.
#' @param CVfolds Control for cross-validation if not supplied in \code{data}.
#' @param parallel \code{boolean} parallelize cross-validation process?
#' Parallelisation is over cross-validation folds.
#' @param pckgs if parallel is TRUE then  specify packages required for
#' each worker (e.g. c("data.table) if data stored as such)
#' @param cores if parallel is TRUE then number of available cores
#' @param Sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param SortLimits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param save_models_path Path to save models. Model details and file extension pasted onto this string.
#' @details The returned predictive quantiles are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' @return Quantile forecasts in a \code{MultiQR} object.
#' @keywords Quantile Regression
#' @export
MQR_qreg_mboost <- function(data,
                            formula,
                            quantiles=c(0.25,0.5,0.75),
                            CVfolds=NULL,
                            ...,
                            bc_mstop=100,
                            bc_nu=0.1,
                            w=rep(1,nrow(data)),
                            parallel = F,
                            cores = NULL,
                            pckgs = NULL,
                            Sort=T,
                            SortLimits=NULL,
                            save_models_path=NULL){
  
  warning("function depreciated and may be removed in future updates. Update to ProbCast::mqr_qreg_mboost()")
  
  
  if(parallel){
    
    
    do.call(mqr_qreg_mboost,c(list(data = data, 
                                   formula = formula,
                                   quantiles=quantiles,
                                   cv_folds=CVfolds,
                                   bc_mstop=bc_mstop,
                                   bc_nu=bc_nu,
                                   w=w,
                                   cores = cores,
                                   pckgs = pckgs,
                                   sort = Sort,
                                   sort_limits= SortLimits,
                                   save_models_path=save_models_path),
                                  ...))
    
    
  } else{
    
    
    do.call(mqr_qreg_mboost,c(list(data = data, 
                                   formula = formula,
                                   quantiles=quantiles,
                                   cv_folds=CVfolds,
                                   bc_mstop=bc_mstop,
                                   bc_nu=bc_nu,
                                   w=w,
                                   cores = 1,
                                   pckgs = pckgs,
                                   sort = Sort,
                                   sort_limits= SortLimits,
                                   save_models_path=save_models_path),
                              ...))
    
    
    
  }
  
  
  
  
  
  
}

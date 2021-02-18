#' Multiple Quantile Regression using \code{mboost}
#'
#' This function fits multiple quantile regression models using \code{mboost}, with
#' facilities for cross-validation. \code{mboost} accommodates both generalised
#' additive models, decision trees and other learners. See \code{?mboost} for more
#' details.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}; Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory variables.
#' @param formaula A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right.
#' @param quantiles The quantiles to fit models for.
#' @param cv_folds Control for cross-validation with various options, either:
#' \itemize{
#'  \item the column name of the fold index supplied in data. Observations and inputs 
#'   in the index labelled "Test" will serve as test data and held out in model training.
#'  \item an integer giving the number of cross validation folds to generate. Folds are constructed as block chunks. 
#'  Default behaviour is 5 folds.
#'  \item vector of length==nrow(data) containing character or numeric fold labels.
#'  \item NULL indicates that no cross validation should be performed and the returned model is trained on all \code{data}.
#' }
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
#' @param only_mqr return only the out-of-sample predictions? 
#' @param exclude_train control for exclusion of rows in data for the model training only, with various options, either:
#' \itemize{
#'  \item the column name of the binary/boolean exclude flag supplied in data.
#'  \item a vector of binary/boolean exclusion flags of length nrow(data)
#'  \item NULL indicates no exclusion
#' }
#' @param ... extra hyper-parameters to be passed to \code{mboost()}.
#' e.g. use \code{control = mboost::boost_control()} to specify boosting steps and shrinkage.
#' @details The returned predictive quantiles are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' @return Quantile forecasts in a \code{MultiQR} object.
#' @keywords Quantile Regression
#' @importFrom foreach %dopar%
#' @export
qreg_mboost <- function(data,
                        formula,
                        quantiles = c(0.25,0.5,0.75),
                        cv_folds = NULL,
                        w = rep(1,nrow(data)),
                        cores = 1,
                        pckgs = NULL,
                        sort = T,
                        sort_limits = NULL,
                        save_models_path = NULL,
                        only_mqr = FALSE,
                        exclude_train = NULL,
                        ...){
  
  ## to-do
  
  if(nrow(data)!=length(w)){
    stop("nrow(data)!=length(w)")
    
  }
  
  if(is.null(cv_folds) & only_mqr){
    stop("no cross validation via cv_folds and return only_mqr is TRUE")
  }
  
  
  output <- list()
  output$call <- match.call()
  
  # set-up cv folds & do checks
  cv_labs <- cv_control(data = data,cv_folds = cv_folds)
  output$kfold_index <- cv_labs$idx
  output$model_names <- cv_labs$fold_loop
  
  # exclude points from training? & do checks
  output$exclude_index <- exclude_fun(data = data,exclude_train = exclude_train)
  
  
  # set up parallel workers, defaults to one worker....
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  #set up progress bar
  iterations <- length(quantiles)
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  gc()
  
  
  
  # fit the models: returns list quantile --> kfold
  output$models <- foreach::foreach(q = quantiles,.packages = c("mboost",pckgs),.options.snow = opts) %dopar% {
    
    temp_qgam <- list()
    
    for(fold in output$model_names){
      
      temp <- mboost::mboost(formula=formula,
                             data=data[output$kfold_index!=fold & output$kfold_index!="Test" & output$exclude_index==0 & !is.na(data[[formula[[2]]]]),],
                             family = mboost::QuantReg(tau=q),
                             weights=w[output$kfold_index!=fold & output$kfold_index!="Test" & output$exclude_index==0 & !is.na(data[[formula[[2]]]])],
                             ...)
      
      if(!is.null(save_models_path)){
        try(save(temp,file = paste0(save_models_path,"_q",100*q,"_fold",fold,".rda")))
      }
      
      temp_qgam[[fold]] <- temp
      
    }
    
    
    return(temp_qgam)
    
    
  }
  
  
  close(pb)
  parallel::stopCluster(cl)
  names(output$models) <- paste0("q",100*quantiles)
  
  
  # re-arrange structure so  kfold --> quantiles
  output$models <- lapply(output$model_names,function(x){
    foldmodel <- lapply(output$models,function(z){
      z[[x]]
    })
    return(foldmodel)
  })
  names(output$models) <- output$model_names
  
  #set new class & default model for prediction
  class(output) <- "qreg_mboost"
  output$default_model <- if("Test"%in%output$model_names){"Test"}else{output$model_names[1]}
  
  
  
  # get cross validation mqr predictions, unless cv_folds = NULL
  output$mqr_pred <- NULL
  if(!is.null(cv_folds)){
    
    #create container for cv output
    output$mqr_pred  <- data.frame(matrix(NA,ncol = length(quantiles), nrow = nrow(data)))
    colnames(output$mqr_pred ) <- paste0("q",100*quantiles)
    
    for(fold in output$model_names){
      
      output$mqr_pred[output$kfold_index==fold,] <- predict(output,
                                                            newdata = data[output$kfold_index==fold,],
                                                            quantiles = NULL,
                                                            model_name = fold,
                                                            sort = sort,
                                                            sort_limits = sort_limits)
      
    }
    
    # class of new cv object
    class(output$mqr_pred) <- c("MultiQR","data.frame")
    
    output$mqr_info <- list(sort = sort,
                            sort_limits = sort_limits)
    
  } else{
    # if no cv delete output$kfold_index...
    output$kfold_index <- NULL
  }
  
  
  
  if(only_mqr){
    
    return(output$mqr_pred)
    
  } else{
    
    return(output)
    
  }
  
  
  
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
  
  warning("function depreciated and may be removed in future updates. Update to ProbCast::qreg_mboost()")
  
  if(is.null(CVfolds)){
    
    fold_cont <- "kfold"
    
  } else{
    
    fold_cont <- CVfolds
    
  }
  
  
  
  if(parallel){
    
    
    do.call(mqr_qreg_mboost,c(list(data = data, 
                                   formula = formula,
                                   quantiles=quantiles,
                                   cv_folds=fold_cont,
                                   control = mboost::boost_control(mstop = bc_mstop,
                                                                   nu = bc_nu),
                                   w=w,
                                   cores = cores,
                                   pckgs = pckgs,
                                   sort = Sort,
                                   sort_limits= SortLimits,
                                   save_models_path=save_models_path,
                                   only_mqr = TRUE),
                                  ...))
    
    
  } else{
    
    
    do.call(mqr_qreg_mboost,c(list(data = data, 
                                   formula = formula,
                                   quantiles=quantiles,
                                   cv_folds=fold_cont,
                                   control = mboost::boost_control(mstop = bc_mstop,
                                                                   nu = bc_nu),
                                   w=w,
                                   cores = 1,
                                   pckgs = pckgs,
                                   sort = Sort,
                                   sort_limits= SortLimits,
                                   save_models_path=save_models_path,
                                   only_mqr = TRUE),
                              ...))
    
    
    
  }
  
  
  
  
  
  
}



#' Predict method for Multiple Quantile Regression using \code{mboost}
#'
#' This function returns multiple quantile predictions as an object of class \code{MultiQR}
#' based on a ProbCast mboost fitted model. S3 Method for for \code{qreg_mboost} objects
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param object 	object of class \code{qreg_mboost} obtained from the function \code{qreg_mboost()}.
#' @param newdata data.frame of observations for which to make predictions
#' @param quantiles the quantiles to predict. Default is all the quantiles present in \code{object}
#' @param model_name the name of the model to be used for prediction.
#' Unless specified, \code{object$default_model} is used.
#' @param sort sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and 
#' lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param ... additional arguments passed on \code{mboost::predict.mboost()}
#' @details this function returns predictive quantiles for each row in \code{newdata},
#' the result is returned as a \code{MultiQR} object 
#' @return Quantile forecasts in a \code{MultiQR} object.
#' @keywords Quantile Regression
#' @export
predict.qreg_mboost <- function(object,
                                newdata = NULL,
                                quantiles = NULL,
                                model_name = NULL,
                                sort = T,
                                sort_limits = NULL,
                                ...) {
  
  
  
  if(class(object)!="qreg_mboost"){stop("object of wrong class, expecting \"qreg_mboost\"")}
  
  
  if(is.null(model_name)){
    model_name <- object$default_model
  } else{ if(sum(model_name%in%object$model_names)!=1){
    
    stop(paste0(model_name," not in names(object)"))
    
  }
  }
  
  
  if(is.null(quantiles)){
    quantiles <- names(object$models[[model_name]])
  }
  
  if(is.numeric(quantiles)){
    quantiles <- paste0("q",quantiles*100)
  }
  
  if(sum(quantiles%in%names(object$models[[model_name]]))!=length(quantiles)){
    stop("specified quantiles not in model")
  }
  
  
  
  pred <- data.frame(sapply(quantiles,function(q){
    
    ### Save out-of-sample predictions
    mboost::predict.mboost(object$models[[model_name]][[q]],,
                           newdata = newdata,
                           ...)
      
    
    
  }))
  
  class(pred) <- c("MultiQR","data.frame")
  
  
  if((sort) & (length(quantiles) != 1)){
    pred <- SortQuantiles(data = pred,Limits = sort_limits)
  }
  
  
  return(pred)
  
  
  
}




### same @export problem as gbm :s


#' Print \code{qreg_mboost} object
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param x object of class \code{qreg_mboost} obtained from the function \code{qreg_mboost()}.
#' @param ... additional arguments; not currently used.
#' @details this function prints details of the call to \code{qreg_mboost()},
#' @keywords Quantile Regression
#' @method print qreg_mboost
#' @export    
print.qreg_mboost <- function(x, ...){
  cat("\n")
  cat("\t Multiple Quantile regression via mboost \n")
  cat("\n")
  cat("Call:\n", deparse(x$call), "\n\n", sep = "")
  
  cat("\n")
  cat("Names of models fitted: ", x$model_names,"\n")
  cat("Default prediction model: ", x$default_model,"\n")
  cat("Quantiles fitted for each model: ", names(x$models[[x$default_model]]),"\n")
  cat("\n")
  invisible(x)
}


### this can be the same as the gbm one

#' Summary \code{qreg_mboost} object
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param x object of class \code{qreg_mboost} obtained from the function \code{qreg_mboost()}.
#' @param ... additional arguments; not currently used.
#' @details this function gives a more detailed description details of the \code{x} object,
#' than \code{print()}
#' @keywords Quantile Regression
#' @method summary qreg_mboost
#' @export
summary.qreg_mboost <- function(x, ...){
  
  summary.qreg_gbm(x)
  
}



#### update mboost model?



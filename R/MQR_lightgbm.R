#' Multiple Quantile Regression Using Gradient Boosted Decision Trees
#'
#' This function fits multiple boosted quantile regression trees 
#' using \code{lightgbm} implementation with facilities for cross-validation.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}; Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}; Rosemary Tawn, \email{rosemary.tawn@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory variables.
#' @param formula A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right. All terms should be column names in \code{data}
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
#' @param num_iterations Number of boosting iterations supplied to \code{lightgbm}, default 100.
#' @param cores the number of available cores. Defaults to one, i.e. no parallelisation
#' @param sort Sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and 
#' lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param save_models_path Path to save models. Model details and file extension pasted onto this string. 
#' Defaults to \code{NULL}, i.e. no model save.
#' @param only_mqr return only the out-of-sample predictions? 
#' @param exclude_train control for exclusion of rows in data for the model training only, with various options, either:
#' \itemize{
#'  \item the column name of the binary/boolean exclude flag supplied in data.
#'  \item a vector of binary/boolean exclusion flags of length nrow(data)
#'  \item NULL indicates no exclusion
#' }
#' This option is useful when out-of-sample predictions are required in rows which need excluded during model training
#' @param ... Additional arguments passed to \code{lightgbm()}.
#' @details The returned predictive quantiles are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' 
#' The returned models are in a named list corresponding to the model for each fold and 
#' and can be extracted for further prediction or evaluation. See \code{predict.qreg_gbm()}.
#' 
#' @return by default a named list containing fitted models as a list of \code{qreg_gbm} objects, 
#' and out-of-sample cross validation  forecasts as an \code{MultiQR} object. The output list depends on \code{cv_folds}.
#' 
#' Alternatively returns only the out-of-sample cross validation forecasts as an \code{MultiQR} 
#' object when \code{only_mqr} is \code{TRUE}
#' @keywords Quantile Regression
#' @export
qreg_lightgbm <- function(data,
                     formula,
                     quantiles = c(0.25,0.5,0.75),
                     cv_folds = NULL,
                     num_iterations=100,
                     cores = 1,
                     sort = TRUE,
                     sort_limits = NULL,
                     save_models_path = NULL,
                     only_mqr = FALSE,
                     exclude_train = NULL,
                     ...){
  
  # to-do
  ## issue/target/fold indexed mqr object?
  ### docu. difference between cv_folds from probcast and cv.folds from gbm
  #### suppress warnings for OOB tree estimates?
  
  
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
  
  
  # get list of features, and the label (response) column, to later make lightgbm dataset.
  features <- labels(terms(formula))
  response <- as.character(formula)[[2]]
  
  if(!all(features %in% names(data))){
    stop("some of the predictors in formula do not have an associated column in data")
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
  

  # fit the models: returns list quantile --> kfold
  output$models <- foreach::foreach(q = quantiles,.packages ="lightgbm",.options.snow = opts) %dopar% {
    temp_lgbm <- list()
    
    for(fold in output$model_names){
      train_data <- data[output$kfold_index!=fold & output$kfold_index!="Test" & output$exclude_index==0 & !is.na(data[[formula[[2]]]]),]
      X <- as.matrix(train_data[, features])
      dataset <- lgb.Dataset(data=X, label=train_data[,response])
      #lgb_params <- list(..., num_iterations=num_iterations, objective="quantile")
      lgb_params <- list(objective="quantile", num_iterations=num_iterations)
      
      lgb_model = lgb.train(params=lgb_params, data=dataset, alpha=q, num_threads=1)
      
      if(!is.null(save_models_path)){
        try(save(lgb_model,file = paste0(save_models_path,"_q",100*q,"_",fold,".rda")))
      }
      
      temp_lgbm[[fold]] <- lgb_model
      
      }
      #if('try-error' %in% class(lgb_model)){
      #  return(paste0("problem fitting model to node ",q))
      #} else{
      #  return(lgb_model)
    return (temp_lgbm)
  
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
  class(output) <- "qreg_lightgbm"
  output$default_model <- if("Test"%in%output$model_names){"Test"}else{output$model_names[1]}
  
  
  # get cross validation mqr predictions, unless cv_folds = NULL
  output$mqr_pred <- NULL
  if(!is.null(cv_folds)){

    #create container for cv output
    output$mqr_pred  <- data.frame(matrix(NA,ncol = length(quantiles), nrow = nrow(data)))
    colnames(output$mqr_pred ) <- paste0("q",100*quantiles)
    
    for(fold in output$model_names){
      
      output$mqr_pred[output$kfold_index==fold,] <- predict(output,
                                                            newdata = as.matrix(data[output$kfold_index==fold,features]),
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






#' Predict method for Multiple Quantile Regression from Gradient Boosted Decision Trees.
#'
#' This function returns multiple quantile predictions as an object of class \code{MultiQR}
#' based on a ProbCast lightgbm fitted model. S3 Method for for \code{qreg_lightgbm} objects
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}; Rosemary Tawn, \email{rosemary.tawn@@strath.ac.uk}
#' @param object 	object of class \code{qreg_lightgbm} obtained from the function \code{qreg_lightgbm()}.
#' @param newdata matrix of observations for which to make predictions; only columns used as model inputs must be included.
#' @param quantiles the quantiles to predict. Default is all the quantiles present in \code{object}
#' @param model_name the name of the model to be used for prediction.
#' Unless specified, \code{object$default_model} is used.
#' @param sort sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and 
#' lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param ... additional arguments; not currently used.
#' @details this function returns predictive quantiles for each row in \code{newdata},
#' the result is returned as a \code{MultiQR} object 
#' @return Quantile forecasts in a \code{MultiQR} object.
#' @keywords Quantile Regression
#' @export
predict.qreg_lgbm <- function(object,
                             newdata = NULL,
                             quantiles = NULL,
                             model_name = NULL,
                             sort = T,
                             sort_limits = NULL,
                             ...) {
  
  
  require(lightgbm)
  if(class(object)!="qreg_lightgbm"){stop("object of wrong class, expecting \"qreg_lightgbm\"")}
  
  
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
    test <- lightgbm::predict.lgb.Booster(object$models[[model_name]][[q]],
                                          data = newdata)
    test <- predict(object$models[[model_name]][[q]],
                                          data = newdata)
 
    
  }))
  
  class(pred) <- c("MultiQR","data.frame")
  
  
  if((sort) & (length(quantiles) != 1)){
    pred <- SortQuantiles(data = pred,Limits = sort_limits)
  }
  
  
  return(pred)
  
  
  
}


### for some reason @export isn't working directly for me in these two functions :s
# https://stackoverflow.com/questions/23724815/roxygen2-issue-with-exporting-print-method
# had to use @method tag as well

#' Print \code{qreg_gbm} object
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param x object of class \code{qreg_gbm} obtained from the function \code{qreg_gbm()}.
#' @param ... additional arguments; not currently used.
#' @details this function prints details of the call to \code{qreg_gbm()},
#' @keywords Quantile Regression
#' @method print qreg_gbm
#' @export    
print.qreg_gbm <- function(x, ...){
  cat("\n")
  cat("\t Multiple Quantile regression via GBM \n")
  cat("\n")
  cat("Call:\n", deparse(x$call), "\n\n", sep = "")

  cat("\n")
  cat("Names of models fitted: ", x$model_names,"\n")
  cat("Default prediction model: ", x$default_model,"\n")
  cat("Quantiles fitted for each model: ", names(x$models[[x$default_model]]),"\n")
  cat("\n")
  invisible(x)
}


#' Summary \code{qreg_gbm} object
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param x object of class \code{qreg_gbm} obtained from the function \code{qreg_gbm()}.
#' @param ... additional arguments; not currently used.
#' @details this function gives a more detailed description details of the \code{x} object,
#' than \code{print()}
#' @keywords Quantile Regression
#' @method summary qreg_gbm
#' @export
summary.qreg_gbm <- function(x, ...){
  
  print(x)
  
  cat("\n")
  
  if(!is.null(x$mqr_pred)){
    cat("out-of-sample kfold cross validation predictions:\n")
    print(data.table::data.table(x$mqr_pred))
    
    cat("\n")
    cat("class of kfold predictions: ", class(x$mqr_pred), "\n")
    cat("sorted quantile predictions?: ", x$mqr_info$sort, "\n")
    if(x$mqr_info$sort){cat("sort limits: ", deparse(x$mqr_info$sort_limits), "\n")}
    
  }
  
}



#### update gbm model?
# https://github.com/dmlc/xgboost/issues/56#issuecomment-53962722
# https://github.com/dmlc/xgboost/issues/3055




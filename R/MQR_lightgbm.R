#' Multiple Quantile Regression Using Gradient Boosted Decision Trees (lightgbm implementation)
#'
#' This function fits multiple boosted quantile regression trees 
#' using \code{lightgbm} with facilities for cross-validation.
#'
#' @author Gordon McFadzean, \email{gordon.mcfadzean@@tneigroup.com}; Rosemary Tawn, \email{rosemary.tawn@@tneigroup.com}
#' @param data A \code{data.frame} containing target and explanatory variables.
#' @param formula A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right. NOTE any manipulation of terms, eg squaring or interactions, within the formula will fail - only individual, linear terms may be specified.
#' @param categoric_features Either a character vector of feature names, or integer vector of indices, for any categoric terms (NULL if not categoric features included).
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
#' @param cores the number of available cores. Defaults to one, i.e. no parallelisation, although in this case the user
#'  must still specify \code{pckgs} if applicable.
#' @param pckgs specify additional packages required for
#' each worker (e.g. c("data.table") if data stored as such).
#' @param sort Sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and
#' lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param only_mqr return only the out-of-sample predictions?
#' @param exclude_train control for exclusion of rows in data for the model training only, with various options, either:
#' \itemize{
#'  \item the column name of the binary/boolean exclude flag supplied in data.
#'  \item a vector of binary/boolean exclusion flags of length nrow(data)
#'  \item NULL indicates no exclusion
#' }
#' This option is useful when out-of-sample predictions are required in rows which need excluded during model training
#' @param lightgbm_params Additional arguments passed to \code{lightgbm()}: objective='quantile' 
#' and the probability level are automatically included so do not need to be specified here.
#' @param ... Additional arguments - not currently used.
#' @details The returned predictive quantiles are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#'
#' The returned models are in a named list corresponding to the model for each fold and
#' and can be extracted for further prediction or evaluation. See \code{predict.qreg_lightgbm()}.
#'
#' @return by default a named list containing fitted models as a list of \code{qreg_lightgbm} objects,
#' and out-of-sample cross validation  forecasts as an \code{MultiQR} object. The output list depends on \code{cv_folds}.
#'
#' Alternatively returns only the out-of-sample cross validation forecasts as an \code{MultiQR}
#' object when \code{only_mqr} is \code{TRUE}
#' @keywords Quantile Regression
#' @importFrom foreach %dopar%
#' @export
qreg_lightgbm <- function(data,
                          formula,
                          categoric_features = NULL,
                          quantiles = c(0.25,0.5,0.75),
                          cv_folds = NULL,
                          cores = 1,
                          pckgs = NULL,
                          sort = TRUE,
                          sort_limits = NULL,
                          only_mqr = FALSE,
                          exclude_train = NULL,
                          lightgbm_params = NULL,
                          ...){

  if(is.null(cv_folds) & only_mqr){
      stop("no cross validation via cv_folds and return only_mqr is TRUE")
  }

  features <- labels(terms(formula))
  response <- as.character(formula)[[2]]

  output <- list()
  output$call <- match.call()
  output$features <- features
  output$response <- response

  # set-up cv folds & do checks
  cv_labs <- cv_control(data = data,cv_folds = cv_folds)
  output$kfold_index <- cv_labs$idx
  output$model_names <- cv_labs$fold_loop

  # exclude points from training? & do checks
  output$exclude_index <- exclude_fun(data = data, exclude_train = exclude_train)

  output$default_model <- if("Test"%in%cv_labs$fold_loop){"Test"}else{cv_labs$fold_loop[1]}
  output$sort <- sort
  output$sort_limits <- sort_limits
  output$lightgbm_params <- lightgbm_params

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
  output$models <- foreach::foreach(q = quantiles, .packages=c("lightgbm", pckgs),.options.snow = opts) %dopar% {
    temp_lgbm <- list()
    for (fold in output$model_names) {
      train_data <- data[output$kfold_index != fold &
                           output$kfold_index != "Test" &
                           output$exclude_index == 0 &
                           !is.na(data[[formula[[2]]]]),]
      X <- as.matrix(train_data %>% select(all_of(features)))
      
      if(is.null(categoric_features)){
        dataset <- lightgbm::lgb.Dataset(data = X, label = train_data[[response]], free_raw_data = FALSE)
      }else{
        dataset <- lightgbm::lgb.Dataset(data = X, label = train_data[[response]], 
                                         categorical_feature=categoric_features, free_raw_data = FALSE)
      }
      
      lgb_params <- c(list(objective = "quantile", alpha = q),
                      lightgbm_params)
      lgb_model <- lightgbm::lightgbm(params = lgb_params, data = dataset)
      temp_lgbm[[fold]] <- lgb_model$save_model_to_string(NULL)

    }
    #if('try-error' %in% class(lgb_model)){
    #  return(paste0("problem fitting model to node ",q))
    #} else{
    #  return(lgb_model)
    return(temp_lgbm)

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


#' Predict method for Multiple Quantile Regression from lightgbm Gradient Boosted Decision Trees.
#'
#' This function returns multiple quantile predictions as an object of class \code{MultiQR}
#' based on a ProbCast gbm fitted model. S3 Method for for \code{qreg_lightgbm} objects
#'
#' @author Gordon McFadzean, \email{gordon.mcfadzean@@tneigroup.com}
#' @param object fitted model of class \code{qreg_lightgbm} obtained from the function \code{qreg_lightgbm()}.
#' @param newdata data.frame of observations for which to make predictions
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
predict.qreg_lightgbm <- function(object,
                             newdata = NULL,
                             quantiles = NULL,
                             model_name = NULL,
                             sort = T,
                             sort_limits = NULL,
                             ...) {



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
      predict(lightgbm::lgb.load(model_str=object$models[[model_name]][[q]]),
              newdata = as.matrix(newdata %>% select(all_of(object$features))))
  }))

  class(pred) <- c("MultiQR","data.frame")


  if((sort) & (length(quantiles) != 1)){
    pred <- SortQuantiles(data = pred,Limits = sort_limits)
  }


  return(pred)



}


#' Model update (retraining) method for Multiple Quantile Regression from lightgbm Gradient Boosted Decision Trees.
#'
#' This function returns updated fitted \code{qreg_lightgbm} models, which are the 
#' supplied previously fitted models updated with new data. S3 Method for \code{qreg_lightgbm} objects
#'
#' @author Gordon McFadzean, \email{gordon.mcfadzean@@tneigroup.com}
#' @param object fitted model of class \code{qreg_lightgbm} obtained from the function \code{qreg_lightgbm()}.
#' @param newdata data.frame of new observations for which to update the model fit
#' @param model_name the name of the model to be used for prediction.
#' Unless specified, \code{object$default_model} is used.
#' @param pckgs specify additional packages required for
#' each worker (e.g. c("data.table") if data stored as such).
#' @param cores the number of available cores. Defaults to one, i.e. no parallelisation, although in this case the user
#'  must still specify \code{pckgs} if applicable. 
#' @param ... additional arguments; not currently used.
#' @details this function updates the fitted lightgbm model with the supplied new data. 
#' NOTE predictions saved within the lightgbm object are not updated within this function.
#' @return a named list containing fitted models as a list of \code{qreg_lightgbm} objects.
#' @keywords Quantile Regression
#' @export
update.qreg_lightgbm <- function(object,
                                 newdata,
                                 model_name=NULL,
                                 pckgs=NULL,
                                 cores=1){


  ## Check class of object
  if(class(object)[1]!="qreg_lightgbm"){stop("object of wrong class, expecting \"qreg_lightgbm\"")}

  ## Use default model unless specified
  if(is.null(model_name)){
    model_name <- object$default_model
  }
  
  quantiles <- as.numeric(sub("q", "", colnames(object$mqr_pred)))/100
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  #set up progress bar
  iterations <- length(quantiles)
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  gc()

  # fit the models: returns list quantile --> kfold
  object$models[[model_name]] <- foreach::foreach(
    q = quantiles,
    .packages =c("lightgbm", pckgs),
    .options.snow = opts) %dopar% {
    X <- as.matrix(newdata[, .SD, .SDcols = object$features])
    dataset <- lightgbm::lgb.Dataset(data = X,
                                     label = newdata[, get(object$response)],
                                     free_raw_data = FALSE)
    lgb_params <- c(list(objective = "quantile", alpha = q),
                    object$lightgbm_params)

    lgb_model <- lightgbm::lightgbm(params = lgb_params,
                                    data = dataset,
                                    init_model = lightgbm::lgb.load(model_str = object$models[[model_name]][[paste0("q",q*100)]]))
    return(lgb_model$save_model_to_string(NULL))
  }


  close(pb)
  parallel::stopCluster(cl)
  names(object$models[[object$default_model]]) <- paste0("q",100*quantiles)
  return(object)
}




#' Model retraining: S3 Generic Method
#'
#' @author Rosemary Tawn, \email{rosemary.tawn@@tneigroup.com}
#' @param object A fitted model object. Currently supported: \code{qreg_lightgbm}.
#' @param ... Additional arguments.
#' @details This is an S3 method, see specific methods \code{\link{retrain_all.qreg_lightgbm}}
#' for details on functionality.
#' @return A fitted model, periodically retrained over the Test set.
#' @export
retrain_all <- function(object, ...) {
  UseMethod("retrain_all", object)
}


#' Method for repeated model retraining throughout the Test set, based on  Multiple 
#' Quantile Regression from lightgbm Gradient Boosted Decision Trees.
#'
#' This function returns updated fitted \code{qreg_lightgbm} models and predictions, which are the 
#' supplied previously fitted models updated throughout the Test set. S3 Method for \code{qreg_lightgbm} objects
#'
#' @author Gordon McFadzean, \email{gordon.mcfadzean@@tneigroup.com}; Rosemary Tawn, \email{rosemary.tawn@@tneigroup.com}
#' @param object fitted model of class \code{qreg_lightgbm} obtained from the function \code{qreg_lightgbm()}.
#' @param data data.frame of model inputs and observations, including a column containing fold labels: 
#' models will be retrained throughout the subset of the data labelled 'Test'.
#' @param retrain_daily_frequency Frequency, in number of days, which the model will be retrained over the Test set.
#' @param issue_datetime_column Name of the column in data containing the forecast issue times. Column of issue times must be POSIXct or Date.
#' @param cv_folds Name of column in data containing the fold labels. Models will 
#' only be retrained for time points beyond the start of the Test fold.
#' @param pckgs specify additional packages required for
#' each worker (e.g. c("data.table") if data stored as such).
#' @param cores the number of available cores. Defaults to one, i.e. no parallelisation, although in this case the user
#'  must still specify \code{pckgs} if applicable. 
#' @param ... additional arguments; not currently used.
#' @details this function periodically retrains the fitted lightgbm model over the Test set and updates the predictions. 
#' NOTE predictions saved within the lightgbm object are not updated within this function.
#' @return a named list containing fitted models as a list of \code{qreg_lightgbm} objects.
#' @keywords Quantile Regression
#' @export
retrain_all.qreg_lightgbm <- function(object,
                                      data,
                                      retrain_daily_frequency,
                                      issue_datetime_column,
                                      cv_folds,
                                      pckgs=NULL,
                                      cores=1,
                                      ...){
  
  if(class(object)!="qreg_lightgbm"){stop("object of wrong class, expecting \"qreg_lightgbm\"")}
  
  # force data to data.table (code is written assuming data is data.table)
  if (!is(data, "data.table")){
    data <- data.table(data)
  }
  
  retrain_limit <- data[get(cv_folds)!='Test', max(get(issue_datetime_column))]
  first_test_issue <- data[get(cv_folds)=='Test', min(get(issue_datetime_column))]

  while(retrain_limit <= data[,max(get(issue_datetime_column))] - lubridate::days(1)){
    retrain_limit <- first_test_issue + lubridate::days(retrain_daily_frequency)

    newdata <- data[get(issue_datetime_column) >= first_test_issue &
                      get(issue_datetime_column) < retrain_limit]
    if(!nrow(data[get(issue_datetime_column) >= retrain_limit])>0){
      break
    }else{
      print("retraining...")
      print(first_test_issue)
      qreg_object <- update.qreg_lightgbm(object=object,
                                          newdata=newdata,
                                          model_name="Test",
                                          pckgs=pckgs,
                                          cores=cores)

      new_preds <- predict(object=object,
                           newdata=data[get(issue_datetime_column) >= retrain_limit],
                           sort=object$sort,
                           sort_limits=object$sort_limits
                           )
      qreg_object$mqr_pred[which(data[[issue_datetime_column]]>= retrain_limit),
                           names(object$mqr_pred)] <- new_preds[, names(object$mqr_pred)]


    }
    first_test_issue <- first_test_issue + lubridate::days(retrain_daily_frequency)
  }
  return(qreg_object)
}


#' Print \code{qreg_lightgbm} object
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param x object of class \code{qreg_lightgbm} obtained from the function \code{qreg_lightgbm()}.
#' @param ... additional arguments; not currently used.
#' @details this function prints details of the call to \code{qreg_lightgbm()},
#' @keywords Quantile Regression
#' @method print qreg_lightgbm
#' @export    
print.qreg_lightgbm <- function(x, ...){
  cat("\n")
  cat("\t Multiple Quantile regression via lightGBM \n")
  cat("\n")
  cat("Call:\n", deparse(x$call), "\n\n", sep = "")
  
  cat("\n")
  cat("Names of models fitted: ", x$model_names,"\n")
  cat("Default prediction model: ", x$default_model,"\n")
  cat("Quantiles fitted for each model: ", names(x$models[[x$default_model]]),"\n")
  cat("\n")
  invisible(x)
}


#' Summary \code{qreg_lightgbm} object
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param x object of class \code{qreg_lightgbm} obtained from the function \code{qreg_lightgbm()}.
#' @param ... additional arguments; not currently used.
#' @details this function gives a more detailed description details of the \code{x} object,
#' than \code{print()}
#' @keywords Quantile Regression
#' @method summary qreg_lightgbm
#' @export
summary.qreg_lightgbm <- function(x, ...){
  
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
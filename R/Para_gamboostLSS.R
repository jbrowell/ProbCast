#' Fit boosted \code{gamlss}-type semi-parametric models
#'
#' @description This function is a wrapper for the function \code{gamboostLSS}, which
#' fits semi-parametric regression models for predictive distributions with up to
#' four parameters (location, scale, shape1, shape2) via gradient boosting.
#'
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory variables.
#' May optionally contain a column called "kfold" with numbered/labeled
#' folds and "Test" for test data.
#' @param formula A formula or named list of formulas 
#' for the location, scale, shape etc..(see \code{?gamboostLSS})
#' @param cv_folds Control for cross-validation if not supplied in \code{data}.
#' @param family A gamboosLSS family object, which is used to define the
#' distribution and the link functions of the various parameters.
#' @param cores the number of available cores. Defaults to one, i.e. no parallelisation, although in this case the user
#'  must still specify \code{pckgs} if applicable.
#' @param pckgs specify additional packages required for
#' each worker (e.g. c("data.table") if data stored as such).
#' @param save_models_path Path to save models. Model details and file extension pasted onto this string. 
#' @param ... additional arguments passed to \code{gamboostLSS()}, 
#' e.g. use control = mboost::boost_control() to specify boosting steps, shrinkage etc.
#' @return A list of \code{gamboostlss} objects. Each list element
#' corresponds to a cross-validation fold and contains a \code{gamlss} model
#' trained on all other folds.
#' @importFrom foreach %dopar%
#' @export
ppd_gamboostlss <- function(data,
                            formula,
                            cv_folds = NULL,
                            family = gamboostLSS::GaussianLSS(),
                            cores = 1,
                            pckgs = NULL,
                            save_models_path = NULL,
                            ...){
  
  ## to-do
  ## -- clear up auto-cv
  ## -- rolling window
  ## -- update S3 methods - class name might clash with source package
  
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
  
  
  # GAMLSS can't handle NAs...
  if(sum(is.na(data))>0){
    warning("NAs in data => data=na.omit(data) passed to gamboostlss().")
  }
  
  
  ## reduce input dataset to strictly necessary inputs
  nms <- c()
  if (is.list(formula)){
    for (i in 1:length(formula)){
      nms <- c(nms,all.names(as.formula(as.character(unname(formula[i])))))
    }
  } else{
    nms <- c(nms,all.names(formula))
  }
  
  nmsind <- which(colnames(data)%in%c(nms,"kfold"))
  tempdata <- na.omit(eval(parse(text=paste0("data[,c(",paste(nmsind,collapse = ","),")]"))))
  
  
  ### set up parallel workers, defaults to one worker....
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  #set up progress bar
  iterations <- length(unique(tempdata$kfold))
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  gc()
  
  
  
  GAMLSSmodelList <- foreach::foreach(fold = unique(tempdata$kfold),.packages = c("gamboostLSS",pckgs),.options.snow = opts) %dopar% {
    
    
    ### works with both data.table and data.frame
    temp <- gamboostLSS::gamboostLSS(data = tempdata[tempdata$kfold!=fold & tempdata$kfold != "Test",],
                                     formula = formula,
                                     families = family,
                                     ...)
    
    
    if(!is.null(save_models_path)){
      try(save(temp,file = paste0(save_models_path,"_fold",fold,".rda")))
    }
    
    
    return(temp)
    
    
  }
  
  close(pb)
  parallel::stopCluster(cl)
  
  
  names(GAMLSSmodelList) <- unique(tempdata$kfold)
  
  # New class to facilitate S3 methods, e.g. PIT.gamboostLSS
  class(GAMLSSmodelList) <- c("gamboostLSS",class(GAMLSSmodelList))
  return(GAMLSSmodelList)
  
  
  
}










#' Fit boosted \code{gamlss}-type semi-parametric models (depreciated)
#'
#' This function is now depreciated and may be removed in future versions of this package.
#' Use \code{ppd_gamboostlss()} instead.

#'
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory variables.
#' May optionally contain a column called "kfold" with numbered/labeled
#' folds and "Test" for test data.
#' @param formula A formula or list of formulas for differences
#' between formulas for location, scale, shape etc..(see \code{?gamboostLSS})
#' @param families A gamboosLSS family object, which is used to define the
#' distribution and the link functions of the various parameters.
#' @param parallel \code{boolean} parallelize cross-validation process?
#' Parallelisation is over cross-validation folds.
#' @param pckgs if parallel is TRUE then  specify packages required for
#' each worker (e.g. c("data.table) if data stored as such)
#' @param ... Additional arguments passed to \code{gamboostLSS()}.
#' @return A list of \code{gamboostlss} objects. Each list element
#' corresponds to a cross-validation fold and contains a \code{gamlss} model
#' trained on all other folds.
#' @export
Para_gamboostLSS <- function(data,
                             formula,
                             families=gamboostLSS::GaussianLSS(),
                             parallel = F,
                             cores = NULL,
                             pckgs = NULL,
                             ...){
  
  
  
  warning("function depreciated and may be removed in future updates. Update to ProbCast::ppd_gamboostlss()")
  
  
  if(parallel){
    
    
    ppd_gamboostlss(data = data,
                    formula = formula,
                    cv_folds = NULL,
                    family = families,
                    cores = cores,
                    pckgs = pckgs,
                    save_models_path = NULL,
                    ...)
    
  } else{
    
    
    ppd_gamboostlss(data = data,
                    formula = formula,
                    cv_folds = NULL,
                    family = families,
                    cores = 1,
                    pckgs = pckgs,
                    save_models_path = NULL,
                    ...)
    
  }
  
}
  
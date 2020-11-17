#' Fit \code{gamlss}-type semi-parametric models
#'
#' @description This function is a wrapper for the function \code{gamlss}, which
#' fits semi-parametric regression models for predictive distributions with up to
#' four parameters (location, scale, shape1, shape2). 
#'
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}; Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory variables.
#' May optionally contain a column called "\code{kfold}" with numbered/labeled folds
#' and "\code{Test}" for test data.
#' @param formula A formula object with the response on the left of an ~ operator,
#' and the terms, separated by + operators, on the right.
#' @param cv_folds Control for cross-validation if not supplied in \code{data}.
#' @param sigma.formula A formula object for fitting a model to the
#' sigma parameter, as in the formula above.
#' @param nu.formula A formula object for fitting a model to the
#' nu parameter, as in the formula above.
#' @param tau.formula A formula object for fitting a model to the
#' tau parameter, as in the formula above.
#' @param family A gamlss.family object, which is used to define the
#' distribution and the link functions of the various parameters.
#' @param cores the number of available cores. Defaults to one, i.e. no parallelisation, although in this case the user
#'  must still specify \code{pckgs} if applicable.
#' @param pckgs specify additional packages required for
#' each worker (e.g. c("data.table") if data stored as such).
#' @param save_models_path Path to save models. Model details and file extension pasted onto this string. 
#' @param ... Additional arguments passed to \code{gamlss()}.
#' @details See \code{?gamlss} for additional details and options.
#' @return A list of \code{gamlss} models with class \code{PPD}. Each list element
#' corresponds to a cross-validation fold and contains a \code{gamlss} model
#' trained on all other folds.
#' @importFrom foreach %dopar%
#' @export
ppd_gamlss <- function(data,
                       formula,
                       cv_folds=NULL,
                       sigma.formula= ~1,
                       nu.formula = ~1,
                       tau.formula = ~1,
                       family = gamlss.dist::NO(),
                       cores = 1,
                       pckgs = NULL,
                       save_models_path=NULL,
                       ... = NULL){
  
  
  ## to-do
  ## -- clear up auto-cv
  ## -- rolling window
  
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
    warning("NAs in data => data=na.omit(data) passed to gamlss().")
  }
  
  nmsind <- which(colnames(data)%in%c(all.names(formula),all.names(sigma.formula),all.names(nu.formula),all.names(tau.formula),"kfold"))
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
  
                         
                         
  GAMLSSmodelList <- foreach::foreach(fold = unique(tempdata$kfold),.packages = c("gamlss",pckgs),.options.snow = opts,.export = "...") %dopar% {
    
    
    
    ### works with both data.table and data.frame
    temp <- gamlss::gamlss(data = tempdata[tempdata$kfold!=fold & tempdata$kfold!="Test",],
                           formula = formula,
                           sigma.formula = sigma.formula,
                           nu.formula = nu.formula,
                           tau.formula = tau.formula,
                           family = family,
                           ...)
    
    if(!is.null(save_models_path)){
      try(save(temp,file = paste0(save_models_path,"_fold",fold,".rda")))
    }
    
    
    return(temp)
    
    
    
  }
  
  close(pb)
  parallel::stopCluster(cl)
  
  
  names(GAMLSSmodelList) <- unique(tempdata$kfold)

  class(GAMLSSmodelList) <- c("PPD",class(GAMLSSmodelList))
  return(GAMLSSmodelList)
    
  
  
}

## - set ellipsis to NuLL because we need to (potentially) export the "method" argument directly in the foreach loop ,.export="..."
## - if not directly exported, like the other functions (e.g. mqr_qreg_gbm()), the mixed()/RS()/CG() functions are not 
##      captured by the foreach automatic exporter because the functions are not defined in the gamlss package...weird :s...
##      they're defined when the model fitting function is called
## - if we don't set the elipsis to NULL, when 'out of the box' settings are used the export throws an error
## - if we directly define it via the the function arguments it also throws an error!




#' Fit \code{gamlss}-type semi-parametric models (depreciated)
#'
#' @description This function is now depreciated and may be removed in future versions of this package.
#' Use \code{ppd_gamlss()} instead.
#'
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}; Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory variables.
#' May optionally contain a column called "\code{kfold}" with numbered/labeled folds
#' and "\code{Test}" for test data.
#' @param formula A formula object with the response on the left of an ~ operator,
#' and the terms, separated by + operators, on the right.
#' @param sigma.formula A formula object for fitting a model to the
#' sigma parameter, as in the formula above.
#' @param nu.formula A formula object for fitting a model to the
#' nu parameter, as in the formula above.
#' @param tau.formula A formula object for fitting a model to the
#' tau parameter, as in the formula above.
#' @param family A gamlss.family object, which is used to define the
#' distribution and the link functions of the various parameters.
#' @param parallel \code{boolean} parallelize model fitting process? Parallelisation is
#' over cross-validation folds.
#' @param pckgs if parallel is TRUE then  specify packages required
#' for each worker (e.g. c("data.table) if data stored as such)
#' @param cores if parallel is TRUE then number of available cores
#' @param ... Additonal arguments passed to \code{gamlss()}.
#' @details See \code{?gamlss} for additional details and options.
#' @return A list of \code{gamlss} models with class \code{PPD}. Each list element
#' corresponds to a cross-validation fold and contains a \code{gamlss} model
#' trained on all other folds.
#' @export
Para_gamlss <- function(data,formula,
                        sigma.formula= ~1,
                        nu.formula = ~1,
                        tau.formula = ~1,
                        family=gamlss.dist::NO(),
                        parallel = F,
                        cores = NULL,
                        pckgs = NULL,
                        ...){
  
  
  warning("function depreciated and may be removed in future updates. Update to ProbCast::ppd_gamlss()")
  
  
  if(parallel){
    
    ppd_gamlss(data = data, 
                formula = formula,
                cv_folds = NULL,
                sigma.formula= sigma.formula,
                nu.formula = nu.formula,
                tau.formula = tau.formula,
                family=family,
                cores = cores,
                pckgs = pckgs,
                save_models_path=NULL,
                ...)
    
    
    
    
  } else{
    
    
    ppd_gamlss(data = data, 
                formula = formula,
                cv_folds = NULL,
                sigma.formula= sigma.formula,
                nu.formula = nu.formula,
                tau.formula = tau.formula,
                family=family,
                cores = 1,
                pckgs = pckgs,
                save_models_path=NULL,
                ...)
    
  }
  
  
  
}

  
  
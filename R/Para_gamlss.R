#' Fit a gamlss paramertirc forecast model
#'
#' @param data A \code{data.frame} containing target and explanatory variables. May optionally contain a collumn called "kfold" with numbered/labeled folds and "Test" for test data.
#' @param formula A formula object with the response on the left of an ~ operator, and the terms, separated by + operators, on the right.
#' @param sigma.formula A formula object for fitting a model to the sigma parameter, as in the formula above.
#' @param nu.formula A formula object for fitting a model to the nu parameter, as in the formula above.
#' @param tau.formula A formula object for fitting a model to the tau parameter, as in the formula above.
#' @param family A gamlss.family object, which is used to define the distribution and the link functions of the various parameters.
#' @param parallel \code{boolean} parallelize model fitting process?
#' @param pckgs if parallel is TRUE then  specify packages required for each worker (e.g. c("data.table) if data stored as such)
#' @param cores if parallel is TRUE then number of available cores
#' @param ... Additonal arguments passed to \code{gamlss()}.
#' @details Details go here...
#' @return A cumulative density function
#' @export
Para_gamlss <- function(data,formula,
                        sigma.formula= ~1,
                        nu.formula = ~1,
                        tau.formula = ~1,
                        family=NO(),
                        parallel = F,
                        cores = NULL,
                        pckgs = NULL,
                        ...){
  
  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    data$kfold<-1
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }
  
  # GAMLSS can't handle NAs...
  if(sum(is.na(data))>0){
    warning("NAs in data => data=na.omit(data) passed to gamlss().")
  }
  
  nmsind <- which(colnames(data)%in%c(all.names(formula),all.names(sigma.formula),all.names(nu.formula),all.names(tau.formula),"kfold"))
  tempdata <- na.omit(eval(parse(text=paste0("data[,c(",paste(nmsind,collapse = ","),")]"))))
  GAMLSSmodelList <- list()
  
  if(parallel){
    
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    iterations <- length(unique(tempdata$kfold))
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    GAMLSSmodelList <- foreach(fold = unique(tempdata$kfold),.packages = c("gamlss",pckgs),.export="...",.options.snow = opts) %dopar% {
      
      ### works with both data.table and data.frame
      temp <- gamlss(data = tempdata[tempdata$kfold!=fold & tempdata$kfold!="Test",],
                     formula = formula,
                     sigma.formula = sigma.formula,
                     nu.formula = nu.formula,
                     tau.formula = tau.formula,
                     family = family,
                     ...)
      
    }
    close(pb)
    stopCluster(cl)
    names(GAMLSSmodelList) <- unique(tempdata$kfold)
    rm(tempdata)
    
    class(GAMLSSmodelList) <- c("PPD",class(GAMLSSmodelList))
    return(GAMLSSmodelList)
    
  } else{
    ### Training Data: k-fold cross-validation/out-of-sample predictions
    
    
    for(fold in unique(tempdata$kfold)){
      print(fold)
      
      temp <- gamlss(data = tempdata[tempdata$kfold!=fold & tempdata$kfold!="Test",],
                     formula = formula,
                     sigma.formula = sigma.formula,
                     nu.formula = nu.formula,
                     tau.formula = tau.formula,
                     family = family,
                     ...)
      
      GAMLSSmodelList[[fold]] <- temp
      
    }
    rm(tempdata)
    class(GAMLSSmodelList) <- c("PPD",class(GAMLSSmodelList))
    return(GAMLSSmodelList)
    
    
  }
  
  
}
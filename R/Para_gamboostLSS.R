#' Fit boosted \code{gamlss}-type semi-paramertirc models
#'
#' @description This function is a wrapper for the function \code{gamboostLSS}, which
#' fits semi-parametric regression models for predictive disributions with up to
#' four parameters (location, scale, shape1, shape2) via gradient boosting.
#'
#' @param data A \code{data.frame} containing target and explanatory variables.
#' May optionally contain a collumn called "kfold" with numbered/labeled
#' folds and "Test" for test data.
#' @param formula A formula or list of formulas for differences
#' between formulas for location, scale, shape etc..(see \code{?gamboostLSS})
#' @param families A gamboosLSS family object, which is used to define the
#' distribution and the link functions of the various parameters.
#' @param parallel \code{boolean} parallelize cross-validation process?
#' Parallelisation is over cross-validatoin folds.
#' @param pckgs if parallel is TRUE then  specify packages required for
#' each worker (e.g. c("data.table) if data stored as such)
#' @param ... Additonal arguments passed to \code{gamboostLSS()}.
#' @return A list of \code{gamboostlss} objects. Each list element
#' corresponsds to a cross-validation fold and contains a \code{gamlss} model
#' trained on all other folds.
#' @export
Para_gamboostLSS <- function(data,formula,families=GaussianLSS(),parallel = F,cores = NULL,pckgs = NULL,...){
  
  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    data$kfold<-1
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }
  
  # GAMLSS can't handle NAs...
  if(sum(is.na(data))>0){
    warning("NAs in data => data=na.omit(data) passed to gamboostlss().")
  }
  
  modelList <- list()
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
  
  if(parallel){
    
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    iterations <- length(unique(tempdata$kfold))
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    modelList <- foreach(fold = unique(tempdata$kfold),.packages = c("gamboostLSS",pckgs),.options.snow = opts) %dopar% {
      
      ### works with both data.table and data.frame
      temp <- gamboostLSS(data = tempdata[tempdata$kfold!=fold & tempdata$kfold != "Test",],
                          formula = formula,
                          families = families,
                          ...)
      
    }
    close(pb)
    stopCluster(cl)
    names(modelList) <- unique(tempdata$kfold)
    rm(tempdata)
    
    return(modelList)
    
  } else{
    ### Training Data: k-fold cross-validation/out-of-sample predictions
    
    
    for(fold in unique(tempdata$kfold)){
      print(fold)
      
      temp <- gamboostLSS(data = tempdata[tempdata$kfold!=fold & tempdata$kfold != "Test",],
                          formula = formula,
                          families = families,
                          ...)
      
      modelList[[fold]] <- temp
      
    }
    rm(tempdata)
    
    return(modelList)
    
    
  }
  
  
}
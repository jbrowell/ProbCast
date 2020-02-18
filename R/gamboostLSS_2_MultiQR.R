#' Convert gamboostLSS objects to MultiQR, or alternatively return predicted parameters of predictive distribution.
#'
#' @param data A \code{data.frame} containing explanatory variables.
#' @param models A Para_gamboostLSS object.
#' @param quantiles Vector of quantiles to be calculated
#' @param params, return distribution parameters instead?
#' @details Details go here...
#' @return A \code{MultiQR} object derived from gamlss predictive distributions. Alternatively, a matrix containing the parameters of the predictive gamboostLSS distributions.
#' @export
gamboostLSS_2_MultiQR <- function(data,models,quantiles=seq(0.05,0.95,by=0.05),params=F){
  
  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    data$kfold<-1
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }
  
  data <- as.data.frame(data)
  # Initialise containers for parameters and quantile forecasts
  parameters <- matrix(1, nrow = nrow(data), ncol =  length(names(models[[1]])))
  colnames(parameters) <- names(models[[1]])
  
  distFamily <- c()
  
  ### current implementation requires all data from gamboost model (don't think you can turn off anyway)
  # models[[fold]]$mu$baselearner$`bbs(SWH)`$get_names() <----possible alternatives.....
  # gamb$mu$baselearner$names(gamb$mu$coef())
  
  
  for(fold in unique(data$kfold)){
    
    tempdata <- data[,which(colnames(data)%in%c(colnames(attributes(models[[fold]])$data),"kfold"))]
    
    # NAs not allowed in newdata. Flags required to record possition.
    gooddata <- rowSums(is.na(tempdata))==0
    
    tempPred <- predict(models[[fold]], newdata = tempdata[data$kfold==fold & gooddata,],type="response")
    
    for(i in 1:(length(tempPred))){
      parameters[data$kfold==fold & gooddata,i] <- tempPred[[i]]
    }
    
    distFamily <- unique(c(distFamily,attributes(attributes(models[[fold]])$families)$name))
    
  }
  
  if(params){
    return(parameters)
  }
  
  
  if(length(distFamily)!=1){stop("length(distFamily)!=1 - Only a single parametric distribution family is allowed.")}
  
  myqfun <- attributes(attributes(models[[fold]])$families)$qfun
  parameters <- data.frame(parameters)
  
  multipleQuantiles <- matrix(NA, nrow = nrow(data), ncol = length(quantiles))
  for (i in 1:length(quantiles)) {
    parameters$p <- quantiles[i]
    multipleQuantiles[, i] <- do.call(myqfun, parameters)
  }
  
  colnames(multipleQuantiles) <- paste0("q",100*quantiles)
  multipleQuantiles <- as.data.frame(multipleQuantiles)
  class(multipleQuantiles) <- c("MultiQR",class(multipleQuantiles))
  
  return(multipleQuantiles)
  
  
}
#' Produce a MultiQR ojbect from a Parametric Predictive Distribution
#' 
#' @description Produce multiple predictive quantiles as a \code{MultiQR} ojbect
#' from a \code{PPD} object. Alternatively
#' return predicted parameters from the \code{PPD} models.
#' 
#' Note that this function may be superseded by an S3 method in future versions of
#' \code{ProbCast}.
#'
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param data A \code{data.frame} containing explanatory variables required by \code{models}.
#' @param models A \code{PPD} object.
#' @param quantiles Vector of quantiles to be included in the returned \code{MultiQR} object.
#'
#' @details Exact quantiles of the (semi-) parametric models (predictive distributions)
#' are calculated and output as a \code{MultiQR} object. Warnings are thrown if the
#' input \code{data} result in \code{NA} estimates of quantiles, likely as a result of inputs
#' being ourside the allowable range for the \code{PPD} model.
#' @return A \code{MultiQR} object derived from gamlss predictive distributions. Alternatively, a matrix condaining the parameters of the predictive  gamlss distributions.
#' @export

PPD_2_MultiQR <- function(data,models,quantiles=seq(0.05,0.95,by=0.05),params=F){
  
  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    data$kfold<-1
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }
  
  # Initialise containers for parameters and quantile forecasts
  parameters <- matrix(1,nrow=nrow(data),ncol=4)
  colnames(parameters) <- c("mu", "sigma", "nu", "tau")
  
  distFamily <- c()
  for(fold in unique(data$kfold)){
    
    tempdata <- data[,which(colnames(data)%in%c(all.names(models[[1]]$mu.formula),
                                                all.names(models[[1]]$sigma.formula),
                                                all.names(models[[1]]$nu.formula),
                                                all.names(models[[1]]$tau.formula)))]
    
    # NAs not allowed in newdata. Flags required to record possition.
    gooddata <- rowSums(is.na(tempdata))==0
    
    tempPred <- predictAll(object = models[[fold]],
                           newdata = tempdata[data$kfold==fold & gooddata,],
                           data = na.omit(tempdata[data$kfold!=fold & data$kfold!="Test",]))
    
    for(i in 1:(length(tempPred)-1)){
      parameters[data$kfold==fold & gooddata,i] <- tempPred[[i]]
    }
    
    distFamily <- unique(c(distFamily,models[[fold]]$family[1]))
    
  }
  
  if(params){
    return(parameters)
  }
  
  ### Calculate Quanties
  if(length(distFamily)!=1){stop("length(distFamily)!=1 - Only a single parametric distribution family is allowed.")}
  
  multipleQuantiles <- matrix(NA,nrow=nrow(data),ncol=length(quantiles))
  
  input <- list()
  if("mu"%in%names(as.list(args(paste0("q",distFamily))))){
    input$mu=parameters[,1]
  }
  if("sigma"%in%names(as.list(args(paste0("q",distFamily))))){
    input$sigma=parameters[,2]
  }
  if("nu"%in%names(as.list(args(paste0("q",distFamily))))){
    input$nu=parameters[,3]
  }
  if("tau"%in%names(as.list(args(paste0("q",distFamily))))){
    input$tau=parameters[,4]
  }
  
  
  for(i in 1:length(quantiles)){
    input$p <- quantiles[i]
    try1 <- try(multipleQuantiles[,i] <- do.call(paste0("q",distFamily),input))
    if(class(try1)=="try-error"){
      for(j in 1:nrow(multipleQuantiles)){
        try(multipleQuantiles[j,i] <- do.call(paste0("q",distFamily),list(mu=parameters[j,1],
                                                                          sigma=parameters[j,2],
                                                                          nu=parameters[j,3],
                                                                          tau=parameters[j,4],
                                                                          p = quantiles[i])))
        warning("NAs in MultiQR")
      }
    }
  }
  
  
  colnames(multipleQuantiles) <- paste0("q",100*quantiles)
  multipleQuantiles <- as.data.frame(multipleQuantiles)
  class(multipleQuantiles) <- c("MultiQR",class(multipleQuantiles))
  
  return(multipleQuantiles)
  
  
}

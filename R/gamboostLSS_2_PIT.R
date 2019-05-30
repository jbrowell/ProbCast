#' Probability integral transform for Para_gamboostLSS objects
#'
#' This function produces a fan plot of a MultiQR object.
#' @param models A Para_gamboostLSS object.
#' @param data Input data corresponding to \code{qrdata}.
#' @param dist_fun cumulative distribution function corresponging to families specified in gamboostLSS model (see example).
#' @param response_name name of response variable in \code{data} object.
#' @details Details go here...
#' @return The probability integral transform of \code{data} through the predictive distribution defined by a list of gamboostLSS objects.
#' @export
gamboostLSS_2_PIT <- function(models,data,dist_fun,response_name,...){
  
  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    if(length(models)!=1){stop("kfold inconsistent with ppd.")}
    data$kfold<-names(models)
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }
  
  data <- as.data.frame(data)
  
  distFamily <- c()
  
  for(fold in unique(data$kfold)){
    distFamily <- unique(c(distFamily,attributes(attributes(models[[fold]])$families)$name))
  }
  
  if(length(distFamily)!=1){stop("length(distFamily)!=1 - Only a single parametric distribution family is allowed.")}
  
  parameters <- gamboostLSS_2_MultiQR(data=data,models=models,params=T)
  
  
  tempdata <- data[,which(colnames(data)%in%c(colnames(attributes(models[[fold]])$data)))]
  gooddata <- rowSums(is.na(tempdata))==0
  
  input <- list(q=data[[response_name]][gooddata])
  
  if("mu"%in%names(as.list(args(dist_fun)))){
    input$mu=parameters[gooddata,1]
  }
  if("sigma"%in%names(as.list(args(dist_fun)))){
    input$sigma=parameters[gooddata,2]
  }
  if("nu"%in%names(as.list(args(dist_fun)))){
    input$nu=parameters[gooddata,3]
  }
  if("tau"%in%names(as.list(args(dist_fun)))){
    input$tau=parameters[gooddata,4]
  }
  
  X <- rep(NA,nrow(data))
  X[gooddata] <- do.call(dist_fun,input)
  
  
  return(X)
}
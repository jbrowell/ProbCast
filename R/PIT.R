#' Probability integral transform: S3 Generic Method
#'
#' @param distdata An object defining cumulative distributions functions. Currently supported: \code{MultiQR}.
#' @param ... Additional arguments.
#' @export
PIT <- function(distdata,...) {
  UseMethod("PIT",distdata)
}

#' Probability integral transform: S3 Method for for MultiQR
#'
#' This function produces a fan plot of a MultiQR object.
#' @param qrdata A \code{MultiQR} object.
#' @param obs A vector of observations corresponding to \code{qrdata}.
#' @param tails A list of arguments passed to \code{condCDF} defining the tails of the CDF
#' @param ... Additional arguments passed to \code{condCDF}.
#' @details Details go here...
#' @return The probability integral transform of \code{obs} through the predictive distribution defined by \code{qrdata} and interpolation scheme in \code{contCDF}.
#' @export
PIT.MultiQR <- function(qrdata,obs,tails,inverse=FALSE,...){
  
  # if(length(obs)!=nrow(qrdata)){stop("length(obs)!=nrow(qrdata)")}
  
  if(tails$method=="exponential" & is.null(tails$thicknessPL) & is.null(tails$thicknessPR)){
    
    if(!("q50"%in%colnames(qrdata))){stop("q50 required for exponential tails.")}
    
    thickness <- rep(NA,tails$nBins)
    targetquants <- stats::quantile(tails$targetvar,probs = seq(0, 1, 1/tails$nBins),na.rm=T)
    for(i in 1:tails$nBins){
      thickness[i] <- mean(tails$targetvar[which(qrdata$q50>=targetquants[i] & qrdata$q50<targetquants[i+1])],na.rm = T)
    }
    tails$thickparamFunc <- stepfun(seq(0,1,length.out = tails$nBins+1),y=c(0,thickness,1))
    
  }
  
  if(tails$method=="ppd_GUtails" | tails$method=="dyn_exponential"){
    print(paste0(tails$method," method valid for [0,1] target variable scale at the moment..."))
  }
  
  
  if (inverse){
    
    X <- matrix(NA,nrow(qrdata),ncol = ncol(obs))
    for(i in 1:nrow(qrdata)){
      if(is.na(qrdata[i,1])){X[i,] <- NA}else{
        if(tails$method=="ppd_GUtails"){
          X[i,] <- contCDF(quantiles = qrdata[i,],tails = append(tails[-which(names(tails)%in%c("lt","rt"))],list(rt=tails$rt[i,],lt=tails$lt[i,])),inverse = TRUE,...)(as.numeric(obs[i,]))} else{
            X[i,] <- contCDF(quantiles = qrdata[i,],tails = tails,inverse = TRUE,...)(as.numeric(obs[i,]))}
      }
    }
  }else{
    
    X<-rep(NA,nrow(qrdata))
    for(i in 1:nrow(qrdata)){
      if(is.na(obs[i]) | is.na(sum(qrdata[i,]))){X[i] <- NA}else{
        if(tails$method=="ppd_GUtails"){
          X[i] <- contCDF(quantiles = qrdata[i,],tails = append(tails[-which(names(tails)%in%c("lt","rt"))],list(rt=tails$rt[i,],lt=tails$lt[i,])),...)(obs[i])} else{
            X[i] <- contCDF(quantiles = qrdata[i,],tails = tails,...)(obs[i])}
      }
    }
    
    # Impose range [0,1] on PIT transformation and issue warning if used.
    if((sum(na.omit(X)>1) + sum(na.omit(X)<0))>0){
      warning("Boundary [0,1] imposed on PIT transformation. Check tails of marginals.")
      X <- ifelse((X <- ifelse(X<0,0,X))>1,1,X)
    }
    
  }
  
  
  return(X)
}

#' Probability integral transform: S3 Method for for PPD
#'
#' This function produces a fan plot of a MultiQR object.
#' @param ppd A \code{PPD} object.
#' @param data Input data corresponding to \code{qrdata}.
#' @param ... Additional arguments passed to \code{condCDF}.
#' @details Details go here...
#' @return The probability integral transform of \code{data} through the predictive distribution defined by \code{ppd}, a list of gamlss objects.
#' @export
PIT.PPD <- function(ppd,data,inverse=FALSE,inv_probs,...){
  
  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    if(length(ppd)!=1){stop("kfold inconsistent with ppd.")}
    data$kfold<-names(ppd)
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }
  
  distFamily <- c()
  for(fold in unique(data$kfold)){
    distFamily <- unique(c(distFamily,ppd[[fold]]$family[1]))
  }
  
  if(length(distFamily)!=1){stop("length(distFamily)!=1 - Only a single parametric distribution family is allowed.")}
  
  parameters <- PPD_2_MultiQR(data=data,models=ppd,params=T)
  
  keepRows <- !is.na(data[[ppd[[1]]$mu.formula[[2]]]])
  
  input <- list(q=data[[ppd[[1]]$mu.formula[[2]]]][keepRows])
  if("mu"%in%names(as.list(args(paste0("q",distFamily))))){
    input$mu=parameters[keepRows,1]
  }
  if("sigma"%in%names(as.list(args(paste0("q",distFamily))))){
    input$sigma=parameters[keepRows,2]
  }
  if("nu"%in%names(as.list(args(paste0("q",distFamily))))){
    input$nu=parameters[keepRows,3]
  }
  if("tau"%in%names(as.list(args(paste0("q",distFamily))))){
    input$tau=parameters[keepRows,4]
  }
  
  if (inverse){
    X <- matrix(NA,nrow(data),ncol = ncol(inv_probs))
    input$q <- NULL
    input$p <- inv_probs[keepRows,]
    X[keepRows,] <- do.call(paste0("q",distFamily),input)
  } else {
    X <- rep(NA,nrow(data))
    X[keepRows] <- do.call(paste0("p",distFamily),input)
  }
  
  return(X)
}
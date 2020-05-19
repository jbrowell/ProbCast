#' Probability integral transform: S3 Generic Method
#'
#' @param distdata An object defining cumulative distributions functions. Currently supported: \code{MultiQR} and \code{PPD}.
#' @param ... Additional arguments.
#' @details This is an S3 method, see spcific methods \code{\link{PIT.MultiQR}}
#' and \code{\link{PIT.PPD}} for details on functionality.
#' @return The Probability integral transform (or its invers) of
#' data through distributions specified by \code{distdata}.
#' @export
PIT <- function(distdata,...) {
  UseMethod("PIT",distdata)
}

#' Probability integral transform: S3 Method for for MultiQR
#'
#' Transforms data (random variables) through their corresponding cumulative
#' distribution function (specified by a \code{MultiQR} object and tail parameters) in 
#' order to produce variables with a standard uniform distribution.
#' @param qrdata A \code{MultiQR} object.
#' @param obs A vector of observations corresponding to the rows of \code{qrdata}.
#' @param tails A list of arguments passed to \code{condCDF} defining the tails of the CDFs.
#' @param inverse A \code{boolean}. If true, the inverse transformation is appiled, i.e. uniform variable
#' to random variable from original distribution.
#' @param ... Additional arguments passed to \code{condCDF}.
#' @details Boundary imposed throwing a warning to ensure that output is in [0,1].
#' @return The probability integral transform of \code{obs} through the
#' predictive distribution defined by \code{qrdata} and
#' interpolation scheme in \code{contCDF} with \code{tails}.
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
  
  if (inverse){
    
    X <- matrix(NA,nrow(qrdata),ncol = ncol(obs))
    for(i in 1:nrow(qrdata)){
      if(is.na(qrdata[i,1])){X[i,] <- NA}else{
        X[i,] <- contCDF(quantiles = qrdata[i,],tails = tails,inverse = TRUE,...)(as.numeric(obs[i,]))}
    }
  }else{
    
    X<-rep(NA,nrow(qrdata))
    for(i in 1:nrow(qrdata)){
      if(is.na(obs[i]) | is.na(sum(qrdata[i,]))){X[i] <- NA}else{
        X[i] <- contCDF(quantiles = qrdata[i,],tails = tails,...)(obs[i])
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
#' Transforms data (random variables) through their corresponding cumulative
#' distribution function (specified by a \code{PPD} object) in 
#' order to produce variables with a standard uniform distribution.
#' @param ppd A \code{PPD} object (a list of \code{gamlss} objects).
#' @param data Input data corresponding to \code{ppd}.
#' @param inverse A \code{boolean}. If true, the inverse transformation is appiled,
#' i.e. uniform variable.
#' @param inv_probs If \code{inverse}, vector of uniform-distributed variable
#' to be transformed.
#' @return The probability integral transform of \code{data} through
#' the predictive distribution defined by \code{PPD}.
#' @export
PIT.PPD <- function(ppd,data,inverse=FALSE,inv_probs=NULL){
  
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
    
    ## Test?: length(inv_probs)==nrow(data)
    
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

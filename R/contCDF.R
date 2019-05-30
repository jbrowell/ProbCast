#' Continuous CFDs
#'
#' This function generats a function that represents a smooth CDF for each row of a MultiQR object.
#' @param quantiles A single-row \code{MultiQR} object.
#' @param kfolds Fold/test label corresponding to \code{quantiles}.
#' @param method Method of interpolation. If \code{method="linear"} linear interpolation is used between quantiles. For spline interpolation, \code{method=list(name=spline,splinemethod)}, where spline method is passed to \code{splinefun}.
#' @param tails Method for tails...
#' @details Details go here...
#' @return A cumulative densift function
#' @export
contCDF <- function(quantiles,kfold=NULL,inverse=F,
                    method=list(name="spline",splinemethod="monoH.FC"),
                    tails=list(method="interpolate",
                               L=0,U=1)){
  ### TESTING
  # quantiles = test1$pred_mqr[500,]
  # inverse=F
  # # tails=list(method="interpolate",
  # #            L=0,U=1)
  # method="linear"
  # tails=list(method="exponential",
  #            L=0,U=1,nBins=10,
  #            DATA=test1,
  #            ntailpoints=5)
  # kfold=NA
  
  
  # rm(quantiles,inverse,tails,method,LnomP,Lquants,Probs,RnomP,Rquants,kfold,i,train,test,thickparamFunc,thickness,thicknessP)
  
  ###
  
  ### Quantiles
  if(nrow(quantiles)!=1){stop("quantiles must be a single-row MultiQR object.")}
  Probs <- as.numeric(gsub(colnames(quantiles),pattern = "q",replacement = ""))/100
  quantiles <- as.numeric(quantiles)
  
  ### Tails
  if(tails$method=="interpolate"){
    LnomP <- 0
    Lquants <- tails$L
    RnomP <- 1
    Rquants <- tails$U
  }else if(tails$method=="interpolate_dtail1"){
    LnomP <- 0
    Lquants <- quantiles[which(Probs==0.5)] + tails$L
    RnomP <- 1
    Rquants <- quantiles[which(Probs==0.5)] + tails$U
  }else if(tails$method=="exponential"){
    
    if(is.null(tails$thickparamFunc)){
      if(!0.5%in%Probs){stop("q50 required for exponential tails.")}
      
      print("single CDF with exponential tails specified")
      
      ##introduce kfold CV into here for defining thickness parameter....
      
      # if(!0.5%in%Probs){stop("q50 required for exponential tails.")}
      # if(is.null(tails$DATA$kfold)){tails$DATA$kfold<-rep(1,nrow(tails$DATA))}
      #
      # ### Calculate Thickness parameters (kfold and test data accounted for)
      # if(!is.na(kfold)){
      #   train <- tails$DATA$kfold==kfold
      #   test <- tails$DATA$kfold!=kfold & !is.na(tails$DATA$kfold)
      # }else{
      #   train <- !is.na(tails$DATA$kfold)
      #   test <- is.na(tails$DATA$kfold)
      # }
      
      thickness <- rep(NA,tails$nBins)
      targetquants <- stats::quantile(tails$targetvar,probs = seq(0, 1, 1/tails$nBins),na.rm=T)
      for(i in 1:tails$nBins){
        thickness[i] <- mean(tails$targetvar[which(tails$preds$q50>=targetquants[i] & tails$preds$q50<targetquants[i+1])],na.rm = T)
      }
      thickparamFunc <- stepfun(seq(0,1,length.out = tails$nBins+1),y=c(0,thickness,1))
      thicknessP <- thickparamFunc(quantiles[which(Probs==0.5)])
    } else {
      thicknessP <- tails$thickparamFunc(quantiles[which(Probs==0.5)])
    }
    
    ### Calculate tails
    if(is.null(tails$ntailpoints)){tails$ntailpoints <- 5}
    
    ### thicknessP has to be less than min(Probs)
    thicknessP <- min(Probs)*10^(-1-thicknessP/tails$U)
    
    if(thicknessP>=min(Probs)){stop("thicknessP has to be less than min(Probs)")}
    
    # Left Tail
    Lquants <- seq(tails$L+min(quantiles)/tails$ntailpoints,by=min(quantiles)/tails$ntailpoints,length.out = tails$ntailpoints-1)
    LnomP <- thicknessP*exp((Lquants/min(quantiles))*log(min(Probs)/thicknessP))
    Lquants <- c(tails$L,Lquants)
    LnomP<-c(0,LnomP)
    # Right Tail
    Rquants <- rev(tails$U-seq((1-max(quantiles))/tails$ntailpoints,by=(1-max(quantiles))/tails$ntailpoints,length.out = tails$ntailpoints-1))
    RnomP <- 1-thicknessP*exp(((1-Rquants)/(1-max(quantiles)))*log((1-max(Probs))/thicknessP))
    Rquants <- c(Rquants,tails$U)
    RnomP<-c(RnomP,1)
    
    # plot(x=c(Lquants,quantiles,Rquants),y=c(LnomP,Probs,RnomP))
    
  }else{stop("Tail specification not recognised.")}
  
  
  
  ### Options: approxfun, splinefun
  if(method=="spline"){
    method <- list(name="spline",splinemethod="monoH.FC")
  }else if(method=="linear"){
    method <- list(name="linear")
  }else{stop("Interpolation method not recognised.")}
  
  if(method$name=="linear"){
    if(!inverse){
      return(approxfun(x=c(Lquants,quantiles,Rquants),y=c(LnomP,Probs,RnomP),yleft = 0,yright = 1))
    }else{
      return(approxfun(x=c(LnomP,Probs,RnomP),y=c(Lquants,quantiles,Rquants),yleft = Lquants[1],yright = tail(Rquants,1)))
    }
  }else if(method$name=="spline"){
    
    if(!inverse){
      return(splinefun(x=c(Lquants,quantiles,Rquants),y=c(LnomP,Probs,RnomP),
                       method=method$splinemethod))
    }else{
      return(splinefun(x=c(LnomP,Probs,RnomP),y=c(Lquants,quantiles,Rquants),
                       method=method$splinemethod))
    }
  }else{stop("Interpolation method not recognised.")}
}

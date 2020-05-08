#' Continuous CFDs
#'
#' This function generats a function that represents a smooth CDF for each row of a MultiQR object.
#' @param quantiles A single-row \code{MultiQR} object.
#' @param kfolds Fold/test label corresponding to \code{quantiles}.
#' @param method Method of interpolation. If \code{method="linear"} linear interpolation is used between quantiles. For spline interpolation, \code{method=list(name=spline,splinemethod)}, where spline method is passed to \code{splinefun}.
#' @param tails Method for tails, this must be defined as a list with other control parameters. Default is linear interpolation of tails when \code{method="interpolate"} with boundaries at \code{L=0} and \code{U=1}. The number of interpolation points for the tail can be defined via \code{ntailpoints}. Alternative methods include:
#' \itemize{
#'   \item Exponential tails can be specified through \code{method="exponential"}, the user can either supply user defined thickness parameters for the tail via \code{thicknessPL} and \code{thicknessPR}, otherwise a symetrical tail thickness can be defined data-driven by specifying: number of bins \code{nbins}, a MQR object \code{preds}, and the target variable, \code{targetvar}.  
#'   \item Dynamic exponential tails can be specified through \code{method="dyn_exponential"}, where the tail depends on the values for the upper and lower quantile of \code{qrdata}, these are only valid for an input variable scale of \code{[0,1]}.
#'   \item The Generalised Pareto Distribution may be chosen with \code{method="gpd"}, and additional arguments \code{scale_r,shape_r,scale_l,shape_l,tail_qs} specifying the right and left scale and shapre parameters and the incraments beyond the final quantiles to evaluate the gpd at, e.g. \code{tail_qs=seq(0.1,1,by=0.1)} to exctend into the tails by +/- 1 in steps of 0.1.
#' }
#' @details Details go here...
#' @return A cumulative density function
#' @export
contCDF <- function(quantiles,kfold=NULL,inverse=F,
                    method=list(name="spline",splinemethod="monoH.FC"),
                    tails=list(method="interpolate",L=0,U=1)){
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
    
    if(!is.null(tails$thicknessPL) & !is.null(tails$thicknessPR)){
      
      
      if(tails$thicknessPL>=min(Probs)){stop("thicknessPL has to be less than min(Probs)")}
      if(tails$thicknessPR>=min(Probs)){stop("thicknessPR has to be less than min(Probs)")}
      if(is.null(tails$ntailpoints)){tails$ntailpoints <- 5}
      
      # Left Tail
      Lquants <- seq(tails$L+min(quantiles)/tails$ntailpoints,by=min(quantiles)/tails$ntailpoints,length.out = tails$ntailpoints-1)
      LnomP <- tails$thicknessPL*exp((Lquants/min(quantiles))*log(min(Probs)/tails$thicknessPL))
      Lquants <- c(tails$L,Lquants)
      LnomP<-c(0,LnomP)
      # Right Tail
      Rquants <- rev(tails$U-seq((1-max(quantiles))/tails$ntailpoints,by=(1-max(quantiles))/tails$ntailpoints,length.out = tails$ntailpoints-1))
      RnomP <- 1-tails$thicknessPR*exp(((1-Rquants)/(1-max(quantiles)))*log((1-max(Probs))/tails$thicknessPR))
      Rquants <- c(Rquants,tails$U)
      RnomP<-c(RnomP,1)
      
    } else{
      if(is.null(tails$thickparamFunc)){
        if(!0.5%in%Probs){stop("q50 required for exponential tails.")}
        
        print("single CDF with exponential tails specified")
        
        ##introduce kfold CV into here for defining thickness parameter....
        
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
      
    }
  }else if(tails$method=="dyn_exponential"){
    
    if(is.null(tails$ntailpoints)){tails$ntailpoints <- 5}
    
    # function only tested for inputs between zero and 1 at the moment, outwith that need to modify to capture minimum bound of tail
    # watch dividing by 0....
    paraf <- function(rho,x,minq){
      if(rho<minq){rho <- minq}
      x*rho*exp((1/x)^(1/(1-rho))*log(minq/rho))
    }
    
    
    # samplebetween 0-1 to get tail shape
    Lquants <- seq(0,1,length.out = tails$ntailpoints)
    # find nominal probabilities
    LnomP <- paraf(min(quantiles),Lquants,min(Probs))
    # remove max value from tail q[(Pr = min(probs))] and normalize shape between available quantile space
    Lquants <- Lquants[1:(length(Lquants)-1)]*min(quantiles)
    # remove max probability from tail (Pr = min(Probs))
    LnomP<-LnomP[1:(length(LnomP)-1)]
    
    # same for R tail
    Rquants <- seq(0,1,length.out = tails$ntailpoints)
    RnomP <- rev(1-paraf(1-max(quantiles),Rquants,minq = 1-max(Probs)))
    Rquants <- rev(((1-Rquants[2:length(Rquants)])*(1-max(quantiles)))+max(quantiles))
    RnomP<-RnomP[2:length(RnomP)]
    
  } else if(tails$method=="gpd"){
    
    Rquants <- tails$tail_qs+rev(quantiles)[1]
    RnomP <- pgpd(q=Rquants,location=rev(quantiles)[1],shape = tails$shape_r,scale=tails$scale_r)
    RnomP[length(RnomP)] <- 1
    RnomP <- RnomP*(1-rev(Probs)[1])+rev(Probs)[1]
    # plot(Rquants,RnomP)
    
    Lquants <- tails$tail_qs
    LnomP <- rev(1-pgpd(q=Lquants,location=0,shape = tails$shape_l,scale=tails$scale_l))
    LnomP[1] <- 0
    LnomP <- LnomP*Probs[1]
    Lquants <- Lquants-max(Lquants)+quantiles[1]
    # plot(Lquants,LnomP)
    
    # plot(c(Lquants,quantiles,Rquants),c(LnomP,Probs,RnomP))
    
    
  }else{stop("Tail specification not recognised.")}
  
  
  
  ### Options: approxfun, splinefun
  if(length(method)==1){
    if(method=="spline"){
      method <- list(name="spline",splinemethod="monoH.FC")
    }else if(method=="linear"){
      method <- list(name="linear")
    }else{stop("Interpolation method not recognised.")}
  }
  
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



pgpd <- function(q,location=0,scale,shape){
  if(shape!=0){
    1-(1+shape*(q-location)/scale)^(-1/shape)
  }else{
    1-exp(-(q-location)/scale)
  }
}


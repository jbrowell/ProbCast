#' Continuous CDF from \code{MultiQR} object 
#'
#' This function generats a smooth, continuous CDF a given row of a \code{MultiQR}
#' object. Interpolation if performed between quantiles and a range of tail models
#' are available for extrapolating beyond beyond the last estimated upper and lower
#' quantile.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}; Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param quantiles A single-row \code{MultiQR} object.
#' @param kfolds Fold/test label corresponding to \code{quantiles}.
#' @param method Method of interpolation. See details.
#' @param tails Definition of tails. See details.
#' @param ... extra arguments to \code{approxfun} or \code{splinefun}.
#' @details Interpolation between quantiles may be linear of via smooth splines:
#' 
#' Linear interpolation: \code{method="linear"} linear interpolation
#' between quantiles.
#' 
#' Spline interpolation: \code{method=list(name=spline,splinemethod=monoH.FC)}, where spline method is
#' passed to \code{splinefun}. \code{splinefun=monoH.FC} is recommended to guarantee monotonically
#' increasing function.
#' 
#' 
#' Several options are available for specifying distribution tails beyond
#' the final upper and lower quantiles:
#' 
#' Linear extrapolation: \code{tails=list(method="extrapolate",L,U)} value set to \code{L} and \code{U}
#' for probability levels 0 and 1, respectively. If \code{method="extrapolate_dtail1"} then
#' tails are extrapolated to the 50th quantile plus (minus) \code{U} (\code{L}).
#' 
#' Exponential tails: \code{tails=list(method="exponential",thicknessPL,thicknessPR,ntailpoints=5)}
#' the user will either supply user defined thickness parameters for the tail
#' via \code{thicknessPL} and \code{thicknessPR}. The number
#' of tail quantiles to be estimated is set by \code{ntailpoints}, which defaults to 5.
#' Alternatively \code{tails=list(method="exponential",thickparamFunc)} where \code{thickparamFunc}
#' is a function that takes the q50 as an input and returns the thickness parameter.
#' 
#' Dynamic exponential tails: \code{tails=list(method="dyn_exponential",ntailpoints=5)}, where the tail shape
#' is conditional on the values for the upper and lower quantile of \code{qrdata}. This method 
#' currently only supports an input variable scale of \code{[0,1]}. The tail shape moves from linear interpolation
#' when the upper/lower quantile is near the boundary for each respective tail, to a conditional exponential shape.
#' 
#' Generalised Pareto Distribution Tails: \code{tails="gpd", scale_r,shape_r,
#' scale_l,shape_l,tail_qs=seq(0.1,2,by=0.1)} with left (_l) and right (_r) scale and shape parameters.
#' Quantiles are calculated at points defined by the upper (lower) quantile plus (minus)
#' \code{tail_qs}.
#' 
#' @return A cumulative density function of the type produced by \code{splinefun} or
#' \code{approxfun}.
#' @export
contCDF <- function(quantiles,kfold=NULL,inverse=F,
                    method=list(name="spline",splinemethod="monoH.FC"),
                    tails=list(method="extrapolate",L=0,U=1),...){
  
  ### Quantiles
  if(nrow(quantiles)!=1){stop("quantiles must be a single-row MultiQR object.")}
  Probs <- as.numeric(gsub(colnames(quantiles),pattern = "q",replacement = ""))/100
  quantiles <- as.numeric(quantiles)
  
  ### Tails
  if(tails$method=="interpolate" | tails$method=="extrapolate"){
    LnomP <- 0
    Lquants <- tails$L
    RnomP <- 1
    Rquants <- tails$U
  }else if(tails$method=="interpolate_dtail1" | tails$method=="extrapolate_dtail1"){
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
      if(!is.null(tails$thickparamFunc)){
        thicknessP <- tails$thickparamFunc(quantiles[which(Probs==0.5)])
        
        ### Calculate tails
        if(is.null(tails$ntailpoints)){tails$ntailpoints <- 5}
        
        ### thicknessP has to be less than min(Probs)
        thicknessP <- min(Probs)*10^(-1-thicknessP/tails$U)
        
        if(thicknessP>=min(Probs)){stop("thicknessP>=min(Probs): thicknessP has to be less than min(Probs)")}
        
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
        
        
      }else{
        stop("Thickness parameters or thickparamFunc must be specified for tail method=\"exponential\"")
      }
    }
  }else if(tails$method=="dyn_exponential"){
    
    if(is.null(tails$ntailpoints)){tails$ntailpoints <- 5}
    if(tails$ntailpoints<5){
      warning("tails$ntailpoints must be at least 5")
      tails$ntailpoints <- 5
    }
    

    # function only tested for inputs between zero and 1 at the moment, outwith that need to modify to capture minimum bound of tail
    paraf_l <- function(lower_q,y,min_p){
      if(lower_q<min_p){lower_q <- min_p}
      y*lower_q*exp(y^(-1/(1-lower_q))*log(min_p/lower_q))
    }

    # samplebetween 0-1 to get tail shape
    Lquants <- seq(0,1,length.out = tails$ntailpoints)
    # remove 0&1
    Lquants <- Lquants[2:(length(Lquants)-1)]
    # find nominal probabilities
    LnomP <- c(0,paraf_l(min(quantiles),Lquants,min(Probs)))
    # normalize shape between available quantile space
    Lquants <- c(0,Lquants*min(quantiles))

    paraf_r <- function(upper_q,y,max_p){
      if(upper_q>max_p){upper_q <- max_p}
      1-y*(1-upper_q)*exp(y^(-1/upper_q)*log((1-max_p)/(1-upper_q)))
    }

    # same for Upp. tail
    Rquants <- seq(0,1,length.out = tails$ntailpoints)
    Rquants <- Rquants[2:(length(Rquants)-1)]

    RnomP <- c(rev(paraf_r(max(quantiles),Rquants,max(Probs))),1)
    Rquants <- c(rev(((1-Rquants)*(1-max(quantiles)))+max(quantiles)),1)
    
    
    
  } else if(tails$method=="gpd"){
    
    ## GPD distribution function
    pgpd <- function(q,location=0,scale,shape){
      if(shape>0){
        1-(1+shape*(q-location)/scale)^(-1/shape)
      }else if(shape<0){
        0 + (q>=location & q<=location-scale/shape)*(1-(1+shape*(q-location)/scale)^(-1/shape))
      }else if(shape==0){
        1-exp(-(q-location)/scale)
      }
    }
    
    if(is.null(tails$tail_qs)){
      tails$tail_qs <- (1:300/100)*abs(diff(range(quantiles)))
    }
    
    Rquants <- tails$tail_qs+rev(quantiles)[1]
    RnomP <- pgpd(q=Rquants,location=rev(quantiles)[1],shape = tails$shape_r,scale=tails$scale_r)
    RnomP[length(RnomP)] <- 1
    RnomP <- RnomP*(1-rev(Probs)[1])+rev(Probs)[1]
    
    Lquants <- tails$tail_qs
    LnomP <- rev(1-pgpd(q=Lquants,location=0,shape = tails$shape_l,scale=tails$scale_l))
    LnomP[1] <- 0
    LnomP <- LnomP*Probs[1]
    Lquants <- Lquants-max(Lquants)+quantiles[1]
    
    # plot(c(Lquants,Rquants),c(LnomP,RnomP))
    # points(c(quantiles),c(Probs),col="red")
    
    
  }else{stop("Tail method not recognised.")}
  
  
  
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
      return(approxfun(x=c(Lquants,quantiles,Rquants),y=c(LnomP,Probs,RnomP),yleft = 0,yright = 1,...))
    }else{
      return(approxfun(x=c(LnomP,Probs,RnomP),y=c(Lquants,quantiles,Rquants),yleft = Lquants[1],yright = tail(Rquants,1),...))
    }
    
  }else if(method$name=="spline"){
    
    if(!inverse){
      return(splinefun(x=c(Lquants,quantiles,Rquants),y=c(LnomP,Probs,RnomP),
                       method=method$splinemethod,...))
    }else{
      return(splinefun(x=c(LnomP,Probs,RnomP),y=c(Lquants,quantiles,Rquants),
                       method=method$splinemethod,...))
    }
    
  }else{stop("Interpolation method not recognised.")}
}




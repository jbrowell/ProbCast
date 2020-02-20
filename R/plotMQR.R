#' Plot a Quantile Forecast
#'
#' This function produces a fan plot of a MultiQR object.
#' @param plotdata A \code{MultiQR} object to be plotted
#' @param targetTimes A vector of forecast target times corresponding to \code{plotdata}.
#' @param quantiles A charater vector. By default the column names from MultiQR object, i.e. all available quantiles. Alternatively, a subset may be specified here.
#' @param q50_line Should the q50 be plotted as a line? This can mis-lead users who may interpret it as the most likely temporal trajectory.
#' @param Legens Location of legend to be produced. E.g. "topleft"
#' @param ... Additional arguments passed to \code{plot()}.
#' @details Details go here...
#' @return A plot of a \code{MQR}.
#' @keywords Quantile Regression, plot
#' @export
plot.MultiQR <- function(plotdata,targetTimes=NULL,quantiles=colnames(plotdata),ylim="auto",q50_line=F,Legend=NULL,...){
  
  # qs <- colnames(plotdata)
  qs <- quantiles
  
  
  if(ylim=="auto"){
    ylim <- range(plotdata)
    ylim <- ylim + rep(diff(ylim),2)*c(-1,1)*0.1
  }
  
  if(!is.null(targetTimes)){
    if(length(targetTimes)!=nrow(plotdata)){stop("length(targetTimes)!=nrow(plotdata)")}
    plotdata$x <- targetTimes
  }else{
    plotdata$x <- 1:nrow(plotdata)
  }
  
  if(!("q50" %in% qs)){stop("MultiQR without q50 not supported")}
  
  
  plot(plotdata$x,plotdata$q50,type="l",ylim=ylim,...)
  
  if(length(qs)>1){
    for(i in 1:floor(length(qs)/2)){
      polygon(c(plotdata$x,rev(plotdata$x)),
              c(plotdata[[qs[i]]],rev(plotdata[[qs[length(qs)+1-i]]])),
              col=rainbow(5*length(qs))[3*length(qs)-i], border=NA)
    }
  }
  
  if(q50_line){lines(plotdata$x,plotdata$q50,col="white")}
  
  if(!is.null(Legend)){
    qs <- as.numeric(gsub("q","",qs))
    legend(Legend,paste0((rev(qs)-qs)[1:floor(length(qs)/2)],"%"),
           pch=15,bty="n",
           col=rainbow(5*length(qs))[3*length(qs)-1:floor(length(qs)/2)],
           ncol=ceiling(length(qs)/10),title = "Prediction Interval")
  }
  
}

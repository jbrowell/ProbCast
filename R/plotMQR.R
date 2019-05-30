#' Plot a Quantile Forecast
#'
#' This function produces a fan plot of a MultiQR object.
#' @param plotdata A \code{MultiQR} object to be plotted
#' @param targetTimes A vector of forecast target times corresponding to \code{plotdata}.
#' @param quantiles A charater vector. By default the column names from MultiQR object, i.e. all available quantiles. Alternatively, a subset may be specified here.
#' @param ... Additional arguments passed to \code{plot()}.
#' @details Details go here...
#' @return A plot of a \code{MQR}.
#' @keywords Quantile Regression, plot
#' @export
plot.MultiQR <- function(plotdata,targetTimes=NULL,quantiles=colnames(plotdata),...){
  
  # qs <- colnames(plotdata)
  qs <- quantiles
  
  if(!is.null(targetTimes)){
    if(length(targetTimes)!=nrow(plotdata)){stop("length(targetTimes)!=nrow(plotdata)")}
    plotdata$x <- targetTimes
  }else{
    plotdata$x <- 1:nrow(plotdata)
  }
  
  if(!("q50" %in% qs)){stop("MultiQR without q50 not supported yet.")}
  
  
  plot(plotdata$x,plotdata$q50,type="l",...)
  
  if(length(qs)>1){
    for(i in 1:floor(length(qs)/2)){
      polygon(c(plotdata$x,rev(plotdata$x)),
              c(plotdata[[qs[i]]],rev(plotdata[[qs[length(qs)+1-i]]])),
              col=rainbow(5*length(qs))[3*length(qs)-i], border=NA)
    }
  }
  
  lines(plotdata$x,plotdata$q50,col="white")
  
}
#' Re-scale MultiQR object
#'
#' This function re-scales MultiQR object.
#' @param mqr A \code{MultiQR} object to be re-scaled
#' @param mult Multiplicative adjustment
#' @param add Additive adjustment
#' @param order Order in which to perform adjustment
#' @details Details go here...
#' @return A \code{MultiQR} object
#' @keywords Quantile Regression, plot
#' @export
rescaleMQR <- function(mqr,mult=1,add=0,Order="ma"){
  
  cl <- class(mqr)
  
  if(Order=="ma"){
    mqr <- mqr*mult+add
  }else if(Order=="am"){
    mqr <- (mqr+add)*mult
  }
  
  class(mqr) <- cl
  
  return(mqr)
  
}
#' Re-scale MultiQR object
#'
#' This function re-scales MultiQR object.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param mqr A \code{MultiQR} object to be re-scaled
#' @param mult Multiplicative adjustment of length 1 or with length equal
#' to the number of rows in mqr.
#' @param add Additive adjustment of length 1 or with
#' length equal to the number of rows in mqr.
#' @param order Order in which to perform adjustment, either
#' \code{"ma"} (multiple-then-add) or vice versa.
#' @details This function re-scales \code{MultiQR} objects. Typical usage
#' is to return forecasts in the domain of the original target variable if modelling
#' and forecasting has been performed after normalisation/standardisation.  
#' @return A \code{MultiQR} object.
#' @keywords Quantile Regression, plot
#' @export
rescaleMQR <- function(mqr,mult=1,add=0,Order="ma"){
  
  if(length(add)!=1 & length(add)!=nrow(mqr)){stop("Lenght of add and mult must be 1 or the same length as there are rows of mqr.")}
  
  cl <- class(mqr)
  
  if(Order=="ma"){
    mqr <- mqr*mult+add
  }else if(Order=="am"){
    mqr <- (mqr+add)*mult
  }
  
  class(mqr) <- cl
  
  return(mqr)
  
}

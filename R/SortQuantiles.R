#' Sort multiple quantiles to remove quantile crossing
#'
#' This function sorts quantiles so that q1<q2<q3<...
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param data A \code{MultiQR, data.frame} containing multiple quantiles
#' @param Limits A \code{list} of the upper (\code{U}) and lower (\code{L}) limits to apply
#' to quantiles if the output should be bounded. E.g. \code{list(U=0.999,L=0.001)}/
#' @details The method use in the produciton of multiple quantiles may not gaurantee
#' that quantiles are (nonstrictly) increasing, which is necseary if combining quantiles
#' to form a cumulative distribution function. Similarly, the quantiles may not
#' respect the boundaries of the target variable. This function re-orders quantiles and 
#' constrains them to the user-specified limits.
#' @return A \code{MultiQR, data.frame} object with (nonstrictly) increasing quanitles.
#' @keywords Quantile Regression
#' @export
SortQuantiles <- function(data,Limits=NULL){

  ### Check cols are in correct order
  if(is.unsorted(as.numeric(gsub(colnames(data),pattern = "q",replacement = "")))){
    stop("Columns are not sorted. Check format.")
  }

  temp <- as.matrix(data)
  temp <- t(apply(temp,1,sort))

  if(!is.null(Limits)){
    temp[temp>Limits$U] <- Limits$U
    temp[temp<Limits$L] <- Limits$L
  }

  temp <- data.frame(temp)
  colnames(temp) <- colnames(data)
  class(temp) <- c("MultiQR","data.frame")

  return(temp)

}

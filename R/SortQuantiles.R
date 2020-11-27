#' Sort multiple quantiles to remove quantile crossing
#'
#' This function sorts quantiles so that q1<q2<q3<...
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param data A \code{MultiQR} containing multiple quantiles
#' @param Limits A \code{list} of the upper (\code{U}) and lower (\code{L}) limits to apply
#' to quantiles if the output should be bounded. E.g. \code{list(U=0.999,L=0.001)}/
#' @details The method use in the production of multiple quantiles may not guarantee
#' that quantiles are (non-strictly) increasing, which is necessary if combining quantiles
#' to form a cumulative distribution function. Similarly, the quantiles may not
#' respect the boundaries of the target variable. This function re-orders quantiles and 
#' constrains them to the user-specified limits.
#' @return A \code{MultiQR, data.frame} object with (non-strictly) increasing quanitles.
#' @keywords Quantile Regression
#' @export
SortQuantiles <- function(data,Limits=NULL){

  ### Check cols are in correct order
  if(is.unsorted(as.numeric(gsub(colnames(data),pattern = "q",replacement = "")))){
    stop("Columns are not sorted. Check format.")
  }
  cl <- class(data)
  if(cl[1]!="MultiQR"){stop("class(data)[1] does not equal \"MultiQR\"")}

  temp <- as.matrix(data)
  
  # Need to exclude rows which are all NAs
  all_nas <- which(rowSums(is.na(temp))==ncol(temp))
  
  if(sum(all_nas)!=0){
    
  temp[-all_nas,] <- t(apply(temp[-all_nas,],1,sort))
  
  } else{
    
    temp <- t(apply(temp,1,sort))
    
  }

  if(!is.null(Limits)){
    temp[temp>Limits$U] <- Limits$U
    temp[temp<Limits$L] <- Limits$L
  }

  temp <- data.frame(temp)
  colnames(temp) <- colnames(data)
  class(temp) <- cl

  return(temp)

}

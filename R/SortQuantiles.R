#' Sort multiple quantiles stored in a data.frame/MultiQR object
#'
#' This function sorts quantiles so that q1<q2<q3<... is satisfied.
#' @param data A \code{MultiQR, data.frame} containing quantile forecasts
#' @param Limits List Upper and lower limits to apply to quantiles. E.g. list(U=0.999,L=0.001)
#' @details Details go here...
#' @return A \code{MultiQR, data.frame} object with ordered quanitles.
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

#' Multivariate Gaussian covariance/correlation matrices
#'
#' This function produces a list of Multivariate Gaussian spatial/tempral/spatiotemporal covarance/correlation matrices from PIT transformed variables
#' @param u_data A dataframe of uniform distributed variables.
#' @param kfold A vector of kfold identifiers.
#' @param cov_cor specify either covariance or correlation
#' @param scale scale gaussian variable inputs to have unit variance for covariance matrix? Defaults to FALSE.
#' @details Details go here...
#' @return A list of covariance/correlation matrices corresponding to kfold ids
#' @export
covcor_matrix <- function(u_data,kfold=NULL,cov_cor="covariance",use="pairwise.complete.obs",scale = F ,...){
  
  ### change to autospecify spatial/spatiotemporal from long format marginals?
  if(is.null(kfold)){
    kfold <- rep(1,nrow(u_data))
  }
  
  
  g_data <- as.data.frame(sapply(u_data, qnorm))
  
  
  matList <- list()
  
  if(cov_cor=="covariance"){
    
    for (fold in unique(kfold)) {
      if(length(unique(kfold))>1){
        if(scale){
          temp <- cov(x = sapply(g_data[kfold!=fold & kfold!="Test",],function(x){x/sd(x,na.rm = T)}), use=use,...)
          temp[temp>1] <- 1
        } else{
          temp <- cov(x = g_data[kfold!=fold & kfold!="Test",], use=use,...)
        }
      }else{
        if(scale){
          temp <- cov(x = sapply(g_data,function(x){x/sd(x,na.rm = T)}), use=use,...)
          temp[temp>1] <- 1
        } else{
          temp <- cov(x = g_data, use=use,...)
        }
      }
      matList[[fold]] <- temp
    }
    
  } else {
    
    for (fold in unique(kfold)) {
      if(length(unique(kfold))>1){
        temp <- cor(x = g_data[kfold!=fold & kfold!="Test", ], use=use,...)
      }else{
        temp <- cor(x = g_data, use=use,...)
      }
      matList[[fold]] <- temp
    }
    
  }
  
  return(matList)
}
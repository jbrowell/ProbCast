#' Sharpness/Interval Width for \code{MultiQR} Objects
#' 
#' @description This function calculates the interval width for symetric quantiles in a \code{MultiQR}
#' object. Optionally, results are produced by cross-validation fold or covariate,
#' 95\% confidence intervals are estimated via bootstrap, and results are plotted.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param qrdata \code{MultiQR} object.
#' @param kfolds Optional vector of fold/test labels corresponding to rows of \code{qrdata}.
#' Cannot be used with \code{subsets}.
#' @param bootstrap Calculate this number of boostrap samples
#' to estimate 95\% confdence interval.
#' @param bootstrap Number of boostrap samples used to generate 95\% confidence intervals.
#' @param ... Additional arguments passed to \code{plot()}.
#' @details Missing values in \code{realisations} are handled by \code{na.rm=T} when
#' calculating average exceedence of a given quantile.
#' @return Quantile Score data and, if \code{plot.it=T}, a plot.
#' @export
sharpness <- function(qrdata,kfolds=NULL,bootstrap=NULL,...){
  
  
  qs <- as.numeric(gsub(colnames(qrdata),pattern = "q",replacement = ""))/100
  int <- qs[qs!=0.5]
  if(!all(int[1:(length(int)/2)] - (1-rev(int[length(int)/2 + 1:(length(int)/2)]))<1e-10)){
    stop("Quantiles are not symetric about q50.")
  }
  int <- (rev(int)-int)[1:(length(int)/2)]
  
  if(!is.null(kfolds)){
    total <- "All_cv"} else{
      total <- "All"  
  }
  
  SRP <- data.frame(Interval=int,
                    Width=as.numeric(rep(NA,length(int))),
                    kfold=total,
                    subset=NA,
                    upper=NA,
                    lower=NA)
  for(q in qs[1:length(int)]){
    if(total == "All_cv"){
      SRP$Width[which(qs==q)] <- mean((qrdata[[paste0("q",100*(1-q))]]-qrdata[[paste0("q",100*q)]])[kfolds!="Test"],
                                   na.rm = T)
    
    if(!is.null(bootstrap)){
      bs_data <- rep(NA,bootstrap)
      for(i in 1:bootstrap){
        data_length <- length(qrdata[[paste0("q",100*q)]][kfolds!="Test"])
        i_samp <- sample(1:data_length,size = data_length,replace = T)
        bs_data[i] <- mean((qrdata[[paste0("q",100*(1-q))]]-qrdata[[paste0("q",100*q)]])[kfolds!="Test"][i_samp],
                           na.rm = T)        
      }
      SRP$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
      SRP$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
    }
    
    
    } else{
      SRP$Width[which(int==(1-q*2))] <- mean((qrdata[[paste0("q",100*(1-q))]]-qrdata[[paste0("q",100*q)]]),
                                             na.rm = T)
    
    
    if(!is.null(bootstrap)){
      bs_data <- rep(NA,bootstrap)
      for(i in 1:bootstrap){
        data_length <- length(qrdata[[paste0("q",100*q)]])
        i_samp <- sample(1:data_length,size = data_length,replace = T)
        bs_data[i] <- mean((qrdata[[paste0("q",100*(1-q))]]-qrdata[[paste0("q",100*q)]])[i_samp],
                           na.rm = T)        
      }
      SRP$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
      SRP$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
    }
    
    
    }
  }
  
  # CV folds
  if(!is.null(kfolds)){
    kfolds[is.na(kfolds)]<-"Test"
    
    for(fold in unique(kfolds)){
      
      tempSRP <- data.frame(Interval=int,
                            Width=as.numeric(rep(NA,length(int))),
                            kfold=fold,
                            subset=NA,
                            upper=NA,
                            lower=NA)
      for(q in qs[1:length(int)]){
        tempSRP$Width[which(qs==q)] <- mean((qrdata[[paste0("q",100*(1-q))]]-qrdata[[paste0("q",100*q)]])[kfolds==fold],
                                           na.rm = T)
        
        if(!is.null(bootstrap)){
          bs_data <- rep(NA,bootstrap)
          for(i in 1:bootstrap){
            data_length <- length(qrdata[[paste0("q",100*q)]][kfolds==fold])
            i_samp <- sample(1:data_length,size = data_length,replace = T)
            bs_data[i] <- mean((qrdata[[paste0("q",100*(1-q))]]-qrdata[[paste0("q",100*q)]])[kfolds==fold][i_samp],na.rm = T)        
          }
          tempSRP$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
          tempSRP$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
        }
        
        
      }
      
      SRP <- rbind(SRP,tempSRP)
    }
  }
  
  
  return(SRP)
}

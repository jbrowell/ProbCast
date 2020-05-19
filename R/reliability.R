#' Reliability Diagram for MultiQR
#'
#' @description Calculated empirical exceedence of each quantile and 
#' plots a reliabiltiy diagram for a \code{MultiQR} object.
#' 
#' Optionally, results may be split by cross-validation fold or covariate and/or
#' confidence intervals may be estimated.
#' @param qrdata \code{MultiQR} object.
#' @param realisations Vector of realisations corresponding to
#' rows of \code{qrdata}. Missing data as \code{NA}s accepted.
#' @param kfolds Optional vector of cross-validation fold labels corresponding
#' to rows of \code{qrdata}. Cannot be used with \code{subsets}.
#' @param subsets Optional vector of covariates to bin data by.
#' Breaks between bins are the empirical quantiles of
#' \code{subsets}. Cannot be used with \code{kfolds}.
#' @param breaks Number of quantiles to divide subsets by, results
#' in \code{breaks+1} bins.
#' @param bootstrap Calculate this number of boostrap samples
#' to estimate 95\% confdence interval
#' @param plot.it \code{boolean}. Make a plot?
#' @param ... Additional arguments passed to \code{plot()}.
#' @details Missing values in \code{realisations} are handled by \code{na.rm=T} when
#' calculating average exceedence of a given quantile.
#' @return Reliability data and, if \code{plot.it=TRUE}, a reliability diagram.
#' @export
reliability <- function(qrdata,realisations,kfolds=NULL,subsets=NULL,breaks=4,bootstrap=NULL,plot.it=T,...){
  
  if(nrow(qrdata)!=length(realisations)){stop("nrow(qrdata)!=length(realisations)")}
  if(!is.null(kfolds)){
    if(length(kfolds)!=length(realisations)){stop("!is.null(kfolds) & nrow(kfolds)!=length(realisations)")}
  }
  if(!is.null(subsets)){
    if(length(subsets)!=length(realisations)){stop("!is.null(subsets) & nrow(subsets)!=length(realisations)")}
  }
  if(!is.null(kfolds) & !is.null(subsets)){stop("Only one of subsets and kfolds can !=NULL.")}
  if(breaks<1){stop("breaks must be a positive integer.")}
  
  qs <- as.numeric(gsub(colnames(qrdata),pattern = "q",replacement = ""))/100
  
  
  if(!is.null(kfolds)){
    total <- "All_cv"
  }else{
    total <- "All"  
  }
  
  Rel <- data.frame(Nominal=qs,
                    Empirical=as.numeric(rep(NA,length(qs))),
                    kfold=total,
                    subset=NA,
                    upper=NA,
                    lower=NA)
  for(q in qs){
    if(total == "All_cv"){
      Rel$Empirical[which(qs==q)] <- mean(qrdata[[paste0("q",100*q)]][kfolds!="Test"]>realisations[kfolds!="Test"],na.rm = T)
      
      if(!is.null(bootstrap)){
        bs_data <- rep(NA,bootstrap)
        for(i in 1:bootstrap){
          data_length <- length(qrdata[[paste0("q",100*q)]][kfolds!="Test"])
          i_samp <- sample(1:data_length,size = data_length,replace = T)
          bs_data[i] <- mean(qrdata[[paste0("q",100*q)]][kfolds!="Test"][i_samp]>realisations[kfolds!="Test"][i_samp],na.rm = T)        
        }
        Rel$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
        Rel$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
      }
      
    }else{
      Rel$Empirical[which(qs==q)] <- mean(qrdata[[paste0("q",100*q)]]>realisations,na.rm = T)
      
      if(!is.null(bootstrap)){
        bs_data <- rep(NA,bootstrap)
        for(i in 1:bootstrap){
          data_length <- length(qrdata[[paste0("q",100*q)]])
          i_samp <- sample(1:data_length,size = data_length,replace = T)
          bs_data[i] <- mean(qrdata[[paste0("q",100*q)]][i_samp]>realisations[i_samp],na.rm = T)        
        }
        Rel$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
        Rel$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
      }
    }
  }
  
  if(plot.it){
    plot(c(0,1),c(0,1),lty=2,type="l",asp = 1,
         xlab = "Nominal",
         ylab = "Empirical",...)
  }
  
  # CV folds
  if(!is.null(kfolds)){
    kfolds[is.na(kfolds)]<-"Test"
    
    for(fold in unique(kfolds)){
      
      tempRel <- data.frame(Nominal=qs,
                            Empirical=as.numeric(rep(NA,length(qs))),
                            kfold=fold,
                            subset=NA,
                            upper=NA,
                            lower=NA)
      for(q in qs){
        tempRel$Empirical[which(qs==q)] <- mean(qrdata[[paste0("q",100*q)]][kfolds==fold]>realisations[kfolds==fold],na.rm = T)
      }
      
      Rel <- rbind(Rel,tempRel)
      
      if(plot.it & fold!="Test"){
        lines(tempRel$Nominal,tempRel$Empirical,type="b",col="Grey50",pch=16,cex=0.7)
      }  
    }
  }
  
  ## Subsets
  if(!is.null(subsets)){
    
    break_qs <- quantile(subsets,probs = seq(from = 1/(breaks+1),by=1/(breaks+1),length.out=breaks),na.rm = T)
    break_qs <- c(-Inf,break_qs,Inf)
    for(i in 2:length(break_qs)){
      indexs <- which(subsets>break_qs[i-1] & subsets<=break_qs[i])  
      
      tempRel <- data.frame(Nominal=qs,
                            Empirical=as.numeric(rep(NA,length(qs))),
                            subset=i-1,
                            kfold=NA,
                            upper=NA,
                            lower=NA)
      for(q in qs){
        tempRel$Empirical[which(qs==q)] <- mean(qrdata[[paste0("q",100*q)]][indexs]>realisations[indexs],na.rm = T)
        
        if(!is.null(bootstrap)){
          bs_data <- rep(NA,bootstrap)
          for(j in 1:bootstrap){
            data_length <- length(qrdata[[paste0("q",100*q)]][indexs])
            i_samp <- sample(1:data_length,size = data_length,replace = T)
            bs_data[j] <- mean(qrdata[[paste0("q",100*q)]][indexs][i_samp]>realisations[indexs][i_samp],na.rm = T)
          }
          tempRel$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
          tempRel$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
        }
        
      }
      
      if(plot.it){
        
        if(!is.null(bootstrap)){
          polygon(x = c(tempRel$Nominal,rev(tempRel$Nominal)),
                  y=c(tempRel$upper,rev(tempRel$lower)),
                  col = rainbow(breaks+1,alpha = 0.3)[i-1],border = NA)
        }
        
        lines(tempRel$Nominal,tempRel$Empirical,type="b",col=rainbow(breaks+1)[i-1],pch=16)
      }
      
      Rel <- rbind(Rel,tempRel)
    }
    
  }
  
  if(plot.it){
    grid()
    if(!is.null(subsets)){
      if(breaks==1){
        legend("topleft",c("Ideal",paste0("<=",break_qs[2]),paste0(">",break_qs[2])),
               lty=c(2,rep(1,breaks+1)),
               col=c(1,rainbow(breaks+1)),
               pch=c(NA,rep(16,breaks+1)),bty = "n")
      }else{
        legend("topleft",c("Ideal",paste0(c("<=",paste0(round(break_qs[2:breaks],digits=1)," to "),">"),
                                          round(break_qs[c(2:(breaks+1),breaks+1)],digits=1))),
               lty=c(2,rep(1,breaks+1)),
               col=c(1,rainbow(breaks+1)),
               pch=c(NA,rep(16,breaks+1)),bty = "n")
      }
    }else{
      
      if(!is.null(bootstrap)){
        # lines(Rel$Nominal[Rel$kfold==total],Rel$upper[Rel$kfold==total],type="l",col=4)
        # lines(Rel$Nominal[Rel$kfold==total],Rel$lower[Rel$kfold==total],type="l",col=4)
        polygon(x = c(Rel$Nominal[Rel$kfold==total],rev(Rel$Nominal[Rel$kfold==total])),
                y=c(Rel$upper[Rel$kfold==total],rev(Rel$lower[Rel$kfold==total])),
                col = rgb(0,0,1,alpha = 0.3),border = NA)
        
      }
      
      lines(Rel$Nominal[Rel$kfold==total],Rel$Empirical[Rel$kfold==total],type="b",col=4,pch=16)
      
      if(!is.null(kfolds) & !("Test"%in%kfolds)){
        legend("topleft",c("Ideal","Forecast","CV Folds"),lty=c(2,1,1),col=c(1,4,"Grey50"),pch=c(NA,16,16),bty = "n")
      }else if("Test"%in%kfolds){
        lines(Rel$Nominal[Rel$kfold=="Test"],Rel$Empirical[Rel$kfold=="Test"],type="b",col="red",pch=16,cex=1)
        legend("topleft",c("Ideal","Test",total,"CV Folds"),lty=c(2,1,1,1),col=c(1,"red",4,"Grey50"),pch=c(NA,16,16,16),bty = "n")
      }else{
        legend("topleft",c("Ideal","Forecast"),lty=c(2,1),col=c(1,4),pch=c(NA,16),bty = "n")
      }
    }
  }
  return(Rel)
}

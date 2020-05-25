#' Pinball Loss for MultiQR
#'
#' This function calculates the pinball loss for each quantile in a \code{MultiQR}
#' object. Optionally, results are produced by cross-validation fold or covariate,
#' 95\% confidence intervals are estimated via bootstrap, and results are plotted.
#'  
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param qrdata \code{MultiQR} object.
#' @param realisations Vector of realisations corresponding to rows of \code{qrdata}.
#' \code{NA} accepted where realisations are missing.
#' @param kfolds Optional vector of fold/test labels corresponding to rows of \code{qrdata}.
#' Cannot be used with \code{subsets}.
#' @param subsets Optional vector of covariates to bin data by.
#' Breaks between bins are the empirical quantiles of
#' \code{subsets}. Cannot be used with \code{kfolds}.
#' @param breaks Number of quantiles to divide subsets by, results
#' in \code{breaks+1} bins.
#' @param bootstrap Calculate this number of boostrap samples
#' to estimate 95\% confdence interval.
#' @param plot.it \code{boolean}. Make a plot?
#' @param ... Additional arguments passed to \code{plot()}.
#' @details Missing values in \code{realisations} are handled by \code{na.rm=T} when
#' calculating average exceedence of a given quantile.
#' @return Quantile Score data and, if plot.it=T, a reliability diagram.
#' @export
pinball <- function(qrdata,realisations,kfolds=NULL,plot.it=T,subsets=NULL,breaks=4,bootstrap=NULL,...){
  
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
    total <- "All_cv"} else{
      total <- "All"  
  }
  
  PBL <- data.frame(Quantile=qs,
                    Loss=as.numeric(rep(NA,length(qs))),
                    kfold=total,
                    subset=NA,
                    upper=NA,
                    lower=NA)
  for(q in qs){
    if(total == "All_cv"){
    PBL$Loss[which(qs==q)] <- mean(((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                     (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]))[kfolds!="Test"],
                                   na.rm = T)
    
    if(!is.null(bootstrap)){
      bs_data <- rep(NA,bootstrap)
      for(i in 1:bootstrap){
        data_length <- length(qrdata[[paste0("q",100*q)]][kfolds!="Test"])
        i_samp <- sample(1:data_length,size = data_length,replace = T)
        bs_data[i] <- mean(((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                              (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]))[kfolds!="Test"][i_samp],
                           na.rm = T)        
      }
      PBL$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
      PBL$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
    }
    
    
    } else{
    PBL$Loss[which(qs==q)] <- mean((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                      (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]),na.rm = T)
    
    
    if(!is.null(bootstrap)){
      bs_data <- rep(NA,bootstrap)
      for(i in 1:bootstrap){
        data_length <- length(qrdata[[paste0("q",100*q)]])
        i_samp <- sample(1:data_length,size = data_length,replace = T)
        bs_data[i] <- mean(((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                             (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]))[i_samp],na.rm = T)        
      }
      PBL$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
      PBL$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
    }
    
    
    }
  }
  
  # CV folds
  if(!is.null(kfolds)){
    kfolds[is.na(kfolds)]<-"Test"
    
    for(fold in unique(kfolds)){
      
      tempPBL <- data.frame(Quantile=qs,
                            Loss=as.numeric(rep(NA,length(qs))),
                            kfold=fold,
                            subset=NA,
                            upper=NA,
                            lower=NA)
      for(q in qs){
        tempPBL$Loss[which(qs==q)] <- mean(((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                              (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]))[kfolds==fold],
                                           na.rm = T)
        
        if(!is.null(bootstrap)){
          bs_data <- rep(NA,bootstrap)
          for(i in 1:bootstrap){
            data_length <- length(qrdata[[paste0("q",100*q)]][kfolds==fold])
            i_samp <- sample(1:data_length,size = data_length,replace = T)
            bs_data[i] <- mean(((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                  (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]))[kfolds==fold][i_samp],na.rm = T)        
          }
          tempPBL$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
          tempPBL$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
        }
        
        
      }
      
      PBL <- rbind(PBL,tempPBL)
    }
  }
  
  
  ## Subsets
  if(!is.null(subsets)){
    
    
    if(is.factor(subsets) | is.character(subsets)){
      
      for(i in unique(subsets)){
        indexs <- which(subsets==i)  
        
        
        tempPBL <- data.frame(Quantile=qs,
                              Loss=as.numeric(rep(NA,length(qs))),
                              kfold=NA,
                              subset = i,
                              upper=NA,
                              lower=NA)
        for(q in qs){
          
          tempPBL$Loss[which(qs==q)] <- mean(((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                                (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]))[indexs],
                                             na.rm = T)
          
          if(!is.null(bootstrap)){
            bs_data <- rep(NA,bootstrap)
            for(j in 1:bootstrap){
              data_length <- length(qrdata[[paste0("q",100*q)]][indexs])
              i_samp <- sample(1:data_length,size = data_length,replace = T)
              bs_data[j] <- mean(((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                    (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]))[indexs][i_samp],na.rm = T)
            }
            tempPBL$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
            tempPBL$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
          }
          
        }
        
        
        PBL <- rbind(PBL,tempPBL)
      }
  
      
      
    } else {
      
    break_qs <- quantile(subsets,probs = seq(from = 1/(breaks+1),by=1/(breaks+1),length.out=breaks),na.rm = T)
    break_qs <- c(-Inf,break_qs,Inf)
      

    for(i in 2:length(break_qs)){
      indexs <- which(subsets>break_qs[i-1] & subsets<=break_qs[i])  
      
      
      tempPBL <- data.frame(Quantile=qs,
                            Loss=as.numeric(rep(NA,length(qs))),
                            kfold=NA,
                            subset = i-1,
                            upper=NA,
                            lower=NA)
      for(q in qs){
        
        tempPBL$Loss[which(qs==q)] <- mean(((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                              (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]))[indexs],
             na.rm = T)
        
        if(!is.null(bootstrap)){
          bs_data <- rep(NA,bootstrap)
          for(j in 1:bootstrap){
            data_length <- length(qrdata[[paste0("q",100*q)]][indexs])
            i_samp <- sample(1:data_length,size = data_length,replace = T)
            bs_data[j] <- mean(((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                  (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]))[indexs][i_samp],na.rm = T)
          }
          tempPBL$upper[which(qs==q)] <- quantile(bs_data,probs = 0.975)
          tempPBL$lower[which(qs==q)] <- quantile(bs_data,probs = 0.025)
        }
        
      }
      
      
      PBL <- rbind(PBL,tempPBL)
    }
      
    }
    
  }
  
  
  
  
  if(plot.it){
    
    
    if(!is.null(subsets)){
      
      plot(PBL[which(PBL$subset==1),1:2],type="b",pch=16,
           xlim=c(0,1),
           ylab="Pinball Loss",col="white",...)
      grid()
      
      lvls <- na.omit(unique(PBL$subset))
      
      for(br in seq_along(lvls)){
        
        lines(PBL[which(PBL$subset==lvls[br]),1:2],type="b",pch=16,col=rainbow(length(lvls))[br])
        
        if(!is.null(bootstrap)){
          polygon(x = c(PBL$Quantile[which(PBL$subset==lvls[br])],rev(PBL$Quantile[which(PBL$subset==lvls[br])])),
                  y=c(PBL$upper[which(PBL$subset==lvls[br])],rev(PBL$lower[which(PBL$subset==lvls[br])])),
                  col = rainbow(length(lvls),alpha = .3)[br],border = NA)
          
        }
        
      }
      
      if(is.factor(subsets) | is.character(subsets)){
        
        legend("topleft",lvls,
               lty=c(rep(1,length(lvls))),
               col=c(rainbow(length(lvls))),
               pch=c(rep(16,length(lvls))),bty = "n",cex=.7,ncol = 2)
        
      } else{
      if(breaks==1){
        legend("topleft",c(paste0("<=",break_qs[2]),paste0(">",break_qs[2])),
               lty=c(rep(1,breaks+1)),
               col=c(rainbow(breaks+1)),
               pch=c(rep(16,breaks+1)),bty = "n")
      }else{
        legend("topleft",c(paste0(c("<=",paste0(round(break_qs[2:breaks],digits=1)," to "),">"),
                                          round(break_qs[c(2:(breaks+1),breaks+1)],digits=1))),
               lty=c(rep(1,breaks+1)),
               col=c(rainbow(breaks+1)),
               pch=c(rep(16,breaks+1)),bty = "n")
      }
      }
    
      
    } else{
      
    
    
    plot(PBL[which(PBL$kfold==total),1:2],type="b",pch=16,
         xlim=c(0,1),
         ylab="Pinball Loss",col="blue",...)
    grid()
    
    if(!is.null(bootstrap)){
      polygon(x = c(PBL$Quantile[which(PBL$kfold==total)],rev(PBL$Quantile[which(PBL$kfold==total)])),
              y=c(PBL$upper[which(PBL$kfold==total)],rev(PBL$lower[which(PBL$kfold==total)])),
              col = rgb(0,0,1,alpha = 0.3),border = NA)
      
    }
    
    if(!is.null(kfolds)){
      for(fold in unique(kfolds)){
        if(fold!="Test"){
          
          lines(PBL[which(PBL$kfold==fold),1:2],type="b",pch=16,col="Grey50")
          
          if(!is.null(bootstrap)){
            polygon(x = c(PBL$Quantile[which(PBL$kfold==fold)],rev(PBL$Quantile[which(PBL$kfold==fold)])),
                    y=c(PBL$upper[which(PBL$kfold==fold)],rev(PBL$lower[which(PBL$kfold==fold)])),
                    col = grey(.5,alpha = 0.3),border = NA)
            
          }
          
        } else{
          lines(PBL[which(PBL$kfold==fold),1:2],type="b",pch=16,col="red")
          
          
          if(!is.null(bootstrap)){
            polygon(x = c(PBL$Quantile[which(PBL$kfold==fold)],rev(PBL$Quantile[which(PBL$kfold==fold)])),
                    y=c(PBL$upper[which(PBL$kfold==fold)],rev(PBL$lower[which(PBL$kfold==fold)])),
                    col = rgb(1,0,0,alpha = 0.3),border = NA)
            
          }
          
        }

      }
      lines(PBL[PBL$kfold==total,1:2],type="b",pch=16,col="blue")
    }
    

    
    if(!is.null(kfolds) & !("Test"%in%kfolds)){
      legend("topleft",c(total,"CV Folds"),lty=c(1,1),col=c("blue","Grey50"),pch=c(16,16),bty = "n")
    }else if("Test"%in%kfolds){
      lines(PBL[PBL$kfold=="Test",1:2],type="b",pch=16,col="red")
      legend("topleft",c("Test",total,"CV Folds"),lty=c(1,1,1),col=c("red","blue","Grey50"),pch=c(16,16,16),bty = "n")
    }
  
    
    }
  }
  
  return(PBL)
}

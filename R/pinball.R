#' Pinball Loss for MultiQR
#'
#' This function calculates and plots the pinball loss for a MultiQR object and plots the result
#' @param qrdata \code{MultiQR} object.
#' @param realisations Vector of realisations corresponding to rows of \code{qrdata}. \code{NA} accepted.
#' @param kfolds Optional vector of fold/test labels corresponding to rows of \code{qrdata}.
#' @param plot.it \code{boolean}. Make a plot?
#' @param ... Additional arguments passed to \code{plot()}.
#' @details Details go here...
#' @return Quantile Score data and, if plot.it=T, a reliability diagram.
#' @export
pinball <- function(qrdata,realisations,kfolds=NULL,plot.it=T,subsets=NULL,breaks=4,bootstrap=NULL,...){
  
  if(nrow(qrdata)!=length(realisations)){stop("nrow(qrdata)!=length(realisations)")}
  if(!is.null(kfolds) & nrow(qrdata)!=length(realisations)){stop("!is.null(kfolds) & nrow(qrdata)!=length(realisations)")}
  
  qs <- as.numeric(gsub(colnames(qrdata),pattern = "q",replacement = ""))/100
  
  if(!is.null(kfolds)){
    total <- "All_cv"} else{
      total <- "All"  
  }
  
  PBL <- data.frame(Quantile=qs,
                    Loss=as.numeric(rep(NA,length(qs))),
                    kfold=total)
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
                            kfold=fold)
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
    
    break_qs <- quantile(subsets,probs = seq(from = 1/(breaks+1),by=1/(breaks+1),length.out=breaks),na.rm = T)
    break_qs <- c(-Inf,break_qs,Inf)
    for(i in 2:length(break_qs)){
      indexs <- which(subsets>break_qs[i-1] & subsets<=break_qs[i])  
      
      
      tempPBL <- data.frame(Quantile=qs,
                            Loss=as.numeric(rep(NA,length(qs))),
                            kfold=NA,
                            subset = i-1)
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
  
  
  
  
  if(plot.it){
    plot(PBL[PBL$kfold==total,1:2],type="b",pch=16,
         xlim=c(0,1),
         ylab="Pinball Loss",col="blue",...)
    grid()
    
    if(!is.null(bootstrap)){
      polygon(x = c(PBL$Quantile[PBL$kfold==total],rev(PBL$Quantile[PBL$kfold==total])),
              y=c(PBL$upper[PBL$kfold==total],rev(PBL$lower[PBL$kfold==total])),
              col = rgb(0,0,1,alpha = 0.3),border = NA)
      
    }
    
    if(!is.null(kfolds)){
      for(fold in unique(kfolds)){
        if(fold!="Test"){
          
          lines(PBL[PBL$kfold==fold,1:2],type="b",pch=16,col="Grey50")
          
          if(!is.null(bootstrap)){
            polygon(x = c(PBL$Quantile[PBL$kfold==fold],rev(PBL$Quantile[PBL$kfold==fold])),
                    y=c(PBL$upper[PBL$kfold==fold],rev(PBL$lower[PBL$kfold==fold])),
                    col = grey(.5,alpha = 0.3),border = NA)
            
          }
          
        } else{
          lines(PBL[PBL$kfold==fold,1:2],type="b",pch=16,col="red")
          
          
          if(!is.null(bootstrap)){
            polygon(x = c(PBL$Quantile[PBL$kfold==fold],rev(PBL$Quantile[PBL$kfold==fold])),
                    y=c(PBL$upper[PBL$kfold==fold],rev(PBL$lower[PBL$kfold==fold])),
                    col = rgb(1,0,0,alpha = 0.3),border = NA)
            
          }
          
        }

      }
      lines(PBL[PBL$kfold==total,1:2],type="b",pch=16,col="blue")
    }
    
    if(!is.null(subsets)){
      
      
      for(sub in seq_along(unique(PBL$subset))){

          lines(PBL[PBL$subset==unique(PBL$subset)[sub],1:2],type="b",pch=16,col=rainbow(length(unique(PBL$subset)))[sub])
          
          if(!is.null(bootstrap)){
            polygon(x = c(PBL$Quantile[PBL$subset==unique(PBL$subset)[sub]],rev(PBL$Quantile[PBL$subset==unique(PBL$subset)[sub]])),
                    y=c(PBL$upper[PBL$subset==unique(PBL$subset)[sub]],rev(PBL$lower[PBL$subset==unique(PBL$subset)[sub]])),
                    col = rainbow(length(unique(PBL$subset)),alpha = .3)[sub],border = NA)
            
          }
      
    
      
    }
    
    
    
    
    if(!is.null(kfolds) & !("Test"%in%kfolds)){
      legend("topleft",c(total,"CV Folds"),lty=c(1,1),col=c("blue","Grey50"),pch=c(16,16),bty = "n")
    }else if("Test"%in%kfolds){
      lines(PBL[PBL$kfold=="Test",1:2],type="b",pch=16,col="red")
      legend("topleft",c("Test",total,"CV Folds"),lty=c(1,1,1),col=c("red","blue","Grey50"),pch=c(16,16,16),bty = "n")
    }
  
    
  }
  
  return(PBL)
}

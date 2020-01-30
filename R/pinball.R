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
pinball <- function(qrdata,realisations,kfolds=NULL,plot.it=T,...){
  
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
                                   na.rm = T)} else{
    PBL$Loss[which(qs==q)] <- mean((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                      (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]),na.rm = T)
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
      }
      
      PBL <- rbind(PBL,tempPBL)
    }
  }
  
  if(plot.it){
    plot(PBL[PBL$kfold==total,1:2],type="b",pch=16,
         ylim=c(min(PBL$Loss),max(PBL$Loss)),
         xlim=c(0,1),
         ylab="Pinball Loss",...)
    
    if(!is.null(kfolds)){
      for(fold in unique(kfolds)){
        lines(PBL[PBL$kfold==fold,1:2],type="b",pch=16,col="Grey50")
      }
      lines(PBL[PBL$kfold==total,1:2],type="b",pch=16)
    }
    
  }
  
  return(PBL)
}

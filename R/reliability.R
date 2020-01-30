#' Reliability Diagram for MultiQR
#'
#' This function plots a reliabiltiy diagram for a MultiQR object and returns the plotdata
#' @param qrdata \code{MultiQR} object.
#' @param realisations Vector of realisations corresponding to rows of \code{qrdata}. \code{NA} accepted.
#' @param kfolds Optional vector of fold/test labels corresponding to rows of \code{qrdata}.
#' @param plot.it \code{boolean}. Make a plot?
#' @param ... Additional arguments passed to \code{plot()}.
#' @details Details go here...
#' @return Reliability data and, if plot.it=T, a reliability diagram.
#' @export
reliability <- function(qrdata,realisations,kfolds=NULL,plot.it=T,...){
  
  if(nrow(qrdata)!=length(realisations)){stop("nrow(qrdata)!=length(realisations)")}
  if(!is.null(kfolds) & nrow(qrdata)!=length(realisations)){stop("!is.null(kfolds) & nrow(qrdata)!=length(realisations)")}
  
  qs <- as.numeric(gsub(colnames(qrdata),pattern = "q",replacement = ""))/100
  
  
  if(!is.null(kfolds)){
    total <- "All_cv"} else{
    total <- "All"  
    }
  
  Rel <- data.frame(Nominal=qs,
                    Empirical=as.numeric(rep(NA,length(qs))),
                    kfold=total)
  for(q in qs){
    if(total == "All_cv"){
    Rel$Empirical[which(qs==q)] <- mean(qrdata[[paste0("q",100*q)]][kfolds!="Test"]>realisations[kfolds!="Test"],na.rm = T)
    } else{
      Rel$Empirical[which(qs==q)] <- mean(qrdata[[paste0("q",100*q)]]>realisations,na.rm = T)
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
                            kfold=fold)
      for(q in qs){
        tempRel$Empirical[which(qs==q)] <- mean(qrdata[[paste0("q",100*q)]][kfolds==fold]>realisations[kfolds==fold],na.rm = T)
      }
      if(plot.it){
        if(fold=="Test"){
          # lines(tempRel$Nominal,tempRel$Empirical,type="b",col="red",pch=16,cex=1)
          # Moved later so that is on top...
        }else{
          lines(tempRel$Nominal,tempRel$Empirical,type="b",col="Grey50",pch=16,cex=0.7)
        }
      }
      Rel <- rbind(Rel,tempRel)
    }
  }
  if(plot.it){
    grid()
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
  
  return(Rel)
}
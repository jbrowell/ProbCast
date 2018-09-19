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



#' Multiple Quantile Regression Using GBM
#'
#' This function fits multiple quantile regreesion GBMs with facilities for cross-validation.
#' @param data A \code{data.frame} containing target and explanatory variables. May optionally contain a collumn called "kfold" with numbered/labeled folds and "Test" for test data.
#' @param formaul A \code{formula} object with the response on the left of an ~ operator, and the terms, separated by + operators, on the right
#' @param quantiles The quantiles to fit models for.
#' @param gbm_params List of parameters to be passed to \code{fit.gbm()}.
#' @param CVfolds Control for cross-validation if not supplied in \code{data}.
#' @param perf.plot Plot GBM performance?
#' @param parallel \code{boolean} parallelize cross-validation process?
#' @param pckgs if parallel is TRUE then  specify packages required for each worker (e.g. c("data.table) if data stored as such)
#' @param cores if parallel is TRUE then number of available cores
#' @param Sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param SortLimits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @details Details go here...
#' @return Quantile forecasts in a \code{MultiQR} object.
#' @keywords Quantile Regression
#' @export
MQR_gbm <- function(data,
                    formula,
                    quantiles=c(0.25,0.5,0.75),
                    CVfolds=NULL,
                    gbm_params=list(interaction.depth = 4,
                                    n.trees = 500,
                                    shrinkage = 0.25,
                                    cv.folds = 1,
                                    n.minobsinnode = 15,
                                    bag.fraction = 0.5,
                                    keep.data = F),
                    perf.plot=F,
                    parallel = F,
                    cores = NULL,
                    pckgs = NULL,
                    Sort=T,SortLimits=NULL){
  
  ### Set-up Cross-validation
  TEST<-F # Flag for Training (with CV) AND Test output
  if("kfold" %in% colnames(data)){
    if(!is.null(CVfolds)){warning("Using column \"kfold\" from data. Argument \"CVfolds\" is not used.")}
    
    if("Test" %in% data$kfold){
      TEST<-T
      nkfold <- length(unique(data$kfold))-1
    }else{
      nkfold <- length(unique(data$kfold))
    }
  }else if(is.null(CVfolds)){
    data$kfold <- rep(1,nrow(data))
    nkfold <- 1
  }else{
    data$kfold <- sort(rep(1:CVfolds,length.out=nrow(data)))
    nkfold <- CVfolds
  }
  
  ### Creae Container for output
  predqs <- data.frame(matrix(NA,ncol = length(quantiles), nrow = nrow(data)))
  colnames(predqs) <- paste0("q",100*quantiles)
  
  if(parallel){
    
    for(q in quantiles){ # Loop over quantiles
      print(paste0("q",q*100))
      
      # Calculate the number of cores
      no_cores <- cores
      # Initiate cluster
      cl <- makeCluster(no_cores)
      registerDoSNOW(cl)
      #set up progress bar
      iterations <- length(unique(data$kfold))
      pb <- txtProgressBar(max = iterations, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      # fit each quantiel model change parameters for CV results
      
      
      
      qpred <- foreach(fold = unique(data$kfold),.packages = c("gbm",pckgs),.options.snow = opts,.combine=c) %dopar% {
        
        ### Fit gbm model
        temp_gbm <- do.call(gbm,c(list(formula=formula,data=data[data$kfold!=fold & data$kfold!="Test" & !is.na(data[[formula[[2]]]]),],distribution = list(name="quantile",alpha=q)),gbm_params))
        
        ### Save out-of-sample predictions
        predict.gbm(temp_gbm,
                    newdata = data[data$kfold==fold,],
                    n.trees = gbm.perf(temp_gbm,plot.it = perf.plot))
        
        
      }
      
      close(pb)
      stopCluster(cl)
      
      predqs[[paste0("q",100*q)]] <- qpred
    }
    
    
  } else{
    ### Training Data: k-fold cross-validation/out-of-sample predictions
    for(q in quantiles){ # Loop over quantiles
      for(fold in unique(data$kfold)){# Loop over CV folds and test data
        
        
        ### Fit gbm model
        temp_gbm <- do.call(gbm,c(list(formula=formula,data=data[data$kfold!=fold & data$kfold!="Test" & !is.na(data[[formula[[2]]]]),],distribution = list(name="quantile",alpha=q)),gbm_params))
        
        ### Save out-of-sample predictions
        predqs[[paste0("q",100*q)]][data$kfold==fold] <- predict.gbm(temp_gbm,
                                                                     newdata = data[data$kfold==fold,],
                                                                     n.trees = gbm.perf(temp_gbm,plot.it = perf.plot))
        
        
        ### Store some performance data?
        
      }
    }
  }
  
  class(predqs) <- c("MultiQR","data.frame")
  

  if(Sort){
    predqs <- SortQuantiles(data = predqs,Limits = SortLimits)
  }
  
  
  
  return(predqs)
  
  
}





#' Plot a Quantile Forecast
#'
#' This function produces a fan plot of a MultiQR object.
#' @param plotdata A \code{MultiQR} object to be plotted
#' @param targetTimes A vector of forecast target times corresponding to \code{plotdata}.
#' @param quantiles A charater vector. By default the column names from MultiQR object, i.e. all available quantiles. Alternatively, a subset may be specified here.
#' @param ... Additional arguments passed to \code{plot()}.
#' @details Details go here...
#' @return A plot of a \code{MQR}.
#' @keywords Quantile Regression, plot
#' @export
plot.MultiQR <- function(plotdata,targetTimes=NULL,quantiles=colnames(plotdata),...){

  # qs <- colnames(plotdata)
  qs <- quantiles

  if(!is.null(targetTimes)){
    if(length(targetTimes)!=nrow(plotdata)){stop("length(targetTimes)!=nrow(plotdata)")}
    plotdata$x <- targetTimes
  }else{
    plotdata$x <- 1:nrow(plotdata)
  }

  if(!("q50" %in% qs)){stop("MultiQR without q50 not supported yet.")}


  plot(plotdata$x,plotdata$q50,type="l",...)

  if(length(qs)>1){
    for(i in 1:floor(length(qs)/2)){
      polygon(c(plotdata$x,rev(plotdata$x)),
              c(plotdata[[qs[i]]],rev(plotdata[[qs[length(qs)+1-i]]])),
              col=rainbow(5*length(qs))[3*length(qs)-i], border=NA)
    }
  }

  lines(plotdata$x,plotdata$q50,col="white")

}


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

  Rel <- data.frame(Nominal=qs,
                    Empirical=as.numeric(rep(NA,length(qs))),
                    kfold="All")
  for(q in qs){
    Rel$Empirical[which(qs==q)] <- mean(qrdata[[paste0("q",100*q)]]>realisations,na.rm = T)
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
    lines(Rel$Nominal[Rel$kfold=="All"],Rel$Empirical[Rel$kfold=="All"],type="b",col=4,pch=16)
    if(!is.null(kfolds) & !("Test"%in%kfolds)){
      legend("topleft",c("Ideal","Forecast","CV Folds"),lty=c(2,1,1),col=c(1,4,"Grey50"),pch=c(NA,16,16),bty = "n")
    }else if("Test"%in%kfolds){
      lines(Rel$Nominal[Rel$kfold=="Test"],Rel$Empirical[Rel$kfold=="Test"],type="b",col="red",pch=16,cex=1)
      legend("topleft",c("Ideal","Test","All","CV Folds"),lty=c(2,1,1,1),col=c(1,"red",4,"Grey50"),pch=c(NA,16,16,16),bty = "n")
    }else{
      legend("topleft",c("Ideal","Forecast"),lty=c(2,1),col=c(1,4),pch=c(NA,16),bty = "n")
    }
  }

  return(Rel)
}


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

  PBL <- data.frame(Quantile=qs,
                    Loss=as.numeric(rep(NA,length(qs))),
                    kfold="All")
  for(q in qs){
    PBL$Loss[which(qs==q)] <- mean((realisations-qrdata[[paste0("q",100*q)]])*q*(realisations>=qrdata[[paste0("q",100*q)]])+
                                     (realisations-qrdata[[paste0("q",100*q)]])*(q-1)*(realisations<qrdata[[paste0("q",100*q)]]),
                                   na.rm = T)
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
    plot(PBL[PBL$kfold=="All",1:2],type="b",pch=16,
         ylim=c(min(PBL$Loss),max(PBL$Loss)),
         xlim=c(0,1),
         ylab="Pinball Loss",...)

    if(!is.null(kfolds)){
      for(fold in unique(kfolds)){
        lines(PBL[PBL$kfold==fold,1:2],type="b",pch=16,col="Grey50")
      }
      lines(PBL[PBL$kfold=="All",1:2],type="b",pch=16)
    }

  }

  return(PBL)
}


#' Continuous CFDs
#'
#' This function generats a function that represents a smooth CDF for each row of a MultiQR object.
#' @param quantiles A single-row \code{MultiQR} object.
#' @param kfolds Fold/test label corresponding to \code{quantiles}.
#' @param method Method of interpolation. If \code{method="linear"} linear interpolation is used between quantiles. For spline interpolation, \code{method=list(name=spline,splinemethod)}, where spline method is passed to \code{splinefun}.
#' @param tails Method for tails...
#' @details Details go here...
#' @return A cumulative densift function
#' @export
contCDF <- function(quantiles,kfold=NULL,inverse=F,
                    method=list(name="spline",splinemethod="monoH.FC"),
                    tails=list(method="interpolate",
                               L=0,U=1)){
  ### TESTING
  # quantiles = test1$pred_mqr[500,]
  # inverse=F
  # # tails=list(method="interpolate",
  # #            L=0,U=1)
  # method="linear"
  # tails=list(method="exponential",
  #            L=0,U=1,nBins=10,
  #            DATA=test1,
  #            ntailpoints=5)
  # kfold=NA


  # rm(quantiles,inverse,tails,method,LnomP,Lquants,Probs,RnomP,Rquants,kfold,i,train,test,thickparamFunc,thickness,thicknessP)

  ###

  ### Quantiles
  if(nrow(quantiles)!=1){stop("quantiles must be a single-row MultiQR object.")}
  Probs <- as.numeric(gsub(colnames(quantiles),pattern = "q",replacement = ""))/100
  quantiles <- as.numeric(quantiles)

  ### Tails
  if(tails$method=="interpolate"){
    LnomP <- 0
    Lquants <- tails$L
    RnomP <- 1
    Rquants <- tails$U
  }else if(tails$method=="exponential"){

    if(is.null(tails$thickparamFunc)){
      if(!0.5%in%Probs){stop("q50 required for exponential tails.")}
      
      print("single CDF with exponential tails specified")
      
      ##introduce kfold CV into here for defining thickness parameter....
      
      # if(!0.5%in%Probs){stop("q50 required for exponential tails.")}
      # if(is.null(tails$DATA$kfold)){tails$DATA$kfold<-rep(1,nrow(tails$DATA))}
      # 
      # ### Calculate Thickness parameters (kfold and test data accounted for)
      # if(!is.na(kfold)){
      #   train <- tails$DATA$kfold==kfold
      #   test <- tails$DATA$kfold!=kfold & !is.na(tails$DATA$kfold)
      # }else{
      #   train <- !is.na(tails$DATA$kfold)
      #   test <- is.na(tails$DATA$kfold)
      # }
      
      thickness <- rep(NA,tails$nBins)
      targetquants <- stats::quantile(tails$targetvar,probs = seq(0, 1, 1/tails$nBins))
      for(i in 1:tails$nBins){
        thickness[i] <- mean(tails$targetvar[which(tails$preds$q50>=targetquants[i] & tails$preds$q50<targetquants[i+1])],na.rm = T)
      }
      thickparamFunc <- stepfun(seq(0,1,length.out = tails$nBins+1),y=c(0,thickness,1))
      thicknessP <- thickparamFunc(quantiles[which(Probs==0.5)])
    } else {
    thicknessP <- tails$thickparamFunc(quantiles[which(Probs==0.5)])
    }

    ### Calculate tails
    if(is.null(tails$ntailpoints)){tails$ntailpoints <- 5}

    ### thicknessP has to be less than min(Probs)
    thicknessP <- min(Probs)*10^(-1-thicknessP/tails$U)

    if(thicknessP>=min(Probs)){stop("thicknessP has to be less than min(Probs)")}

    # Left Tail
    Lquants <- seq(tails$L+min(quantiles)/tails$ntailpoints,by=min(quantiles)/tails$ntailpoints,length.out = tails$ntailpoints-1)
    LnomP <- thicknessP*exp((Lquants/min(quantiles))*log(min(Probs)/thicknessP))
    Lquants <- c(tails$L,Lquants)
    LnomP<-c(0,LnomP)
    # Right Tail
    Rquants <- rev(tails$U-seq((1-max(quantiles))/tails$ntailpoints,by=(1-max(quantiles))/tails$ntailpoints,length.out = tails$ntailpoints-1))
    RnomP <- 1-thicknessP*exp(((1-Rquants)/(1-max(quantiles)))*log((1-max(Probs))/thicknessP))
    Rquants <- c(Rquants,tails$U)
    RnomP<-c(RnomP,1)

    # plot(x=c(Lquants,quantiles,Rquants),y=c(LnomP,Probs,RnomP))

  }else{stop("Tail specification not recognised.")}



  ### Options: approxfun, splinefun
  if(method=="spline"){
    method <- list(name="spline",splinemethod="monoH.FC")
  }else if(method=="linear"){
    method <- list(name="linear")
  }else{stop("Interpolation method not recognised.")}

  if(method$name=="linear"){
    if(!inverse){
      return(approxfun(x=c(Lquants,quantiles,Rquants),y=c(LnomP,Probs,RnomP),yleft = 0,yright = 1))
    }else{
      return(approxfun(x=c(LnomP,Probs,RnomP),y=c(Lquants,quantiles,Rquants),yleft = Lquants[1],yright = tail(Rquants,1)))
    }
  }else if(method$name=="spline"){

    if(!inverse){
      return(splinefun(x=c(Lquants,quantiles,Rquants),y=c(LnomP,Probs,RnomP),
                       method=method$splinemethod))
    }else{
      return(splinefun(x=c(LnomP,Probs,RnomP),y=c(Lquants,quantiles,Rquants),
                       method=method$splinemethod))
    }
  }else{stop("Interpolation method not recognised.")}
}




#' Fit a gamlss paramertirc forecast model
#'
#' @param data A \code{data.frame} containing target and explanatory variables. May optionally contain a collumn called "kfold" with numbered/labeled folds and "Test" for test data.
#' @param formula A formula object with the response on the left of an ~ operator, and the terms, separated by + operators, on the right.
#' @param sigma.formula A formula object for fitting a model to the sigma parameter, as in the formula above.
#' @param nu.formula A formula object for fitting a model to the nu parameter, as in the formula above.
#' @param tau.formula A formula object for fitting a model to the tau parameter, as in the formula above.
#' @param family A gamlss.family object, which is used to define the distribution and the link functions of the various parameters.
#' @param ... Additonal arguments passed to \code{gamlss()}.
#' @details Details go here...
#' @return A cumulative densift function
#' @export
Para_gamlss <- function(data,formula,
                        sigma.formula= ~1,
                        nu.formula = ~1,
                        tau.formula = ~1,
                        family=NO(),
                        ...){

  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    data$kfold<-1
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }

  # GAMLSS can't handle NAs...
  if(sum(is.na(data))>0){
    warning("NAs in data => data=na.omit(data) passed to gamlss().")
  }

  GAMLSSmodelList <- list()

  for(fold in unique(data$kfold)){


    temp <- gamlss(formula = formula,
                   sigma.formula = sigma.formula,
                   nu.formula = nu.formula,
                   tau.formula = tau.formula,
                   family = family,
                   data = na.omit(data[data$kfold!=fold & data$kfold!="Test",
                                       which(colnames(data)%in%c(all.names(formula),
                                                                 all.names(sigma.formula),
                                                                 all.names(nu.formula),
                                                                 all.names(tau.formula)))]),
                   ...)

    GAMLSSmodelList[[fold]] <- temp

  }

  class(GAMLSSmodelList) <- c("PPD",class(GAMLSSmodelList))


  return(GAMLSSmodelList)

}


#' Convert PPD to MultiQR, or alternatively return predicted parameters of predictive distribution.
#'
#' @param data A \code{data.frame} containing explanatory variables.
#' @param models An \code{PPD} object.
#' @param quantiles Vector of quantiles to be calculated
#'
#' @details Details go here...
#' @return A \code{MultiQR} object derived from gamlss predictive distributions. Alternatively, a matrix condaining the parameters of the predictive  gamlss distributions.
#' @export

PPD_2_MultiQR <- function(data,models,quantiles=seq(0.05,0.95,by=0.05),params=F){

  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    data$kfold<-1
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }

  # Initialise containers for parameters and quantile forecasts
  parameters <- matrix(1,nrow=nrow(data),ncol=4)
  colnames(parameters) <- c("mu", "sigma", "nu", "tau")

  distFamily <- c()
  for(fold in unique(data$kfold)){

    tempdata <- data[,which(colnames(data)%in%c(all.names(models[[1]]$mu.formula),
                                                all.names(models[[1]]$sigma.formula),
                                                all.names(models[[1]]$nu.formula),
                                                all.names(models[[1]]$tau.formula)))]

    # NAs not allowed in newdata. Flags required to record possition.
    gooddata <- rowSums(is.na(tempdata))==0

    tempPred <- predictAll(object = models[[fold]],
                           newdata = tempdata[data$kfold==fold & gooddata,],
                           data = na.omit(tempdata[data$kfold!=fold & data$kfold!="Test",]))

    for(i in 1:(length(tempPred)-1)){
      parameters[data$kfold==fold & gooddata,i] <- tempPred[[i]]
    }

    distFamily <- unique(c(distFamily,models[[fold]]$family[1]))

  }

  if(params){
    return(parameters)
  }

  ### Calculate Quanties
  if(length(distFamily)!=1){stop("length(distFamily)!=1 - Only a single parametric distribution family is allowed.")}

  multipleQuantiles <- matrix(NA,nrow=nrow(data),ncol=length(quantiles))

  input <- list()
  if("mu"%in%names(as.list(args(paste0("q",distFamily))))){
    input$mu=parameters[,1]
  }
  if("sigma"%in%names(as.list(args(paste0("q",distFamily))))){
    input$sigma=parameters[,2]
  }
  if("nu"%in%names(as.list(args(paste0("q",distFamily))))){
    input$nu=parameters[,3]
  }
  if("tau"%in%names(as.list(args(paste0("q",distFamily))))){
    input$tau=parameters[,4]
  }


  for(i in 1:length(quantiles)){
    input$p <- quantiles[i]
    multipleQuantiles[,i] <- do.call(paste0("q",distFamily),input)

  }


  colnames(multipleQuantiles) <- paste0("q",100*quantiles)
  multipleQuantiles <- as.data.frame(multipleQuantiles)
  class(multipleQuantiles) <- c("MultiQR",class(multipleQuantiles))

  return(multipleQuantiles)


}


#' Probability integral transform: S3 Generic Method
#'
#' @param distdata An object defining cumulative distributions functions. Currently supported: \code{MultiQR}.
#' @param ... Additional arguments.
#' @export
PIT <- function(distdata,...) {
  UseMethod("PIT",distdata)
}

#' Probability integral transform: S3 Method for for MultiQR
#'
#' This function produces a fan plot of a MultiQR object.
#' @param qrdata A \code{MultiQR} object.
#' @param obs A vector of observations corresponding to \code{qrdata}.
#' @param tails A list of arguments passed to \code{condCDF} defining the tails of the CDF
#' @param ... Additional arguments passed to \code{condCDF}.
#' @details Details go here...
#' @return The probability integral transform of \code{obs} through the predictive distribution defined by \code{qrdata} and interpolation scheme in \code{contCDF}.
#' @export
PIT.MultiQR <- function(qrdata,obs,tails,...){

  if(length(obs)!=nrow(qrdata)){stop("length(obs)!=nrow(qrdata)")}
  
  if(tails$method=="exponential"){
    
    if(!("q50"%in%colnames(qrdata))){stop("q50 required for exponential tails.")}
    
    ##introduce kfold CV into here for defining thickness parameter....
    
    # if(!0.5%in%Probs){stop("q50 required for exponential tails.")}
    # if(is.null(tails$DATA$kfold)){tails$DATA$kfold<-rep(1,nrow(tails$DATA))}
    # 
    # ### Calculate Thickness parameters (kfold and test data accounted for)
    # if(!is.na(kfold)){
    #   train <- tails$DATA$kfold==kfold
    #   test <- tails$DATA$kfold!=kfold & !is.na(tails$DATA$kfold)
    # }else{
    #   train <- !is.na(tails$DATA$kfold)
    #   test <- is.na(tails$DATA$kfold)
    # }
    
    thickness <- rep(NA,tails$nBins)
    targetquants <- stats::quantile(tails$targetvar,probs = seq(0, 1, 1/tails$nBins))
    for(i in 1:tails$nBins){
      thickness[i] <- mean(tails$targetvar[which(qrdata$q50>=targetquants[i] & qrdata$q50<targetquants[i+1])],na.rm = T)
    }
    tails$thickparamFunc <- stepfun(seq(0,1,length.out = tails$nBins+1),y=c(0,thickness,1))
    
  }

  X<-rep(NA,nrow(qrdata))
  for(i in 1:nrow(qrdata)){
    if(is.na(obs[i])){X[i] <- NA}else{
    X[i] <- contCDF(quantiles = qrdata[i,],tails = tails,...)(obs[i])}
  }
  return(X)
}

#' Probability integral transform: S3 Method for for PPD
#'
#' This function produces a fan plot of a MultiQR object.
#' @param ppd A \code{PPD} object.
#' @param data Input data corresponding to \code{qrdata}.
#' @param ... Additional arguments passed to \code{condCDF}.
#' @details Details go here...
#' @return The probability integral transform of \code{data} through the predictive distribution defined by \code{ppd}, a list of gamlss objects.
#' @export
PIT.PPD <- function(ppd,data,...){

  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    if(length(ppd)!=1){stop("kfold inconsistent with ppd.")}
    data$kfold<-names(ppd)
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }

  distFamily <- c()
  for(fold in unique(data$kfold)){
    distFamily <- unique(c(distFamily,ppd[[fold]]$family[1]))
  }

  if(length(distFamily)!=1){stop("length(distFamily)!=1 - Only a single parametric distribution family is allowed.")}

  parameters <- PPD_2_MultiQR(data=data,models=ppd,params=T)

  keepRows <- !is.na(data[[ppd[[1]]$mu.formula[[2]]]])

  input <- list(q=data[[ppd[[1]]$mu.formula[[2]]]][keepRows])
  if("mu"%in%names(as.list(args(paste0("q",distFamily))))){
    input$mu=parameters[keepRows,1]
  }
  if("sigma"%in%names(as.list(args(paste0("q",distFamily))))){
    input$sigma=parameters[keepRows,2]
  }
  if("nu"%in%names(as.list(args(paste0("q",distFamily))))){
    input$nu=parameters[keepRows,3]
  }
  if("tau"%in%names(as.list(args(paste0("q",distFamily))))){
    input$tau=parameters[keepRows,4]
  }

  X <- rep(NA,nrow(data))
  X[keepRows] <- do.call(paste0("p",distFamily),input)

  return(X)
}


#' Fit a gamboostlss paramertirc forecast model
#'
#' @param data A \code{data.frame} containing target and explanatory variables. May optionally contain a collumn called "kfold" with numbered/labeled folds and "Test" for test data.
#' @param formula A formula or list of formulas for differences between formulas for location, scale, shape etc..(see \code{gamboostLSS)
#' @param families A gamboosLSS family object, which is used to define the distribution and the link functions of the various parameters.
#' @param parallel \code{boolean} parallelize cross-validation process?
#' @param pckgs if parallel is TRUE then  specify packages required for each worker (e.g. c("data.table) if data stored as such)
#' @param ... Additonal arguments passed to \code{gamboostLSS()}.
#' @return A list of gamboostlss objects corresponsing to kfolds
#' @export
Para_gamboostLSS <- function(data,formula,families=GaussianLSS(),parallel = F,cores = NULL,pckgs = NULL,...){
  
  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    data$kfold<-1
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }
  
  # GAMLSS can't handle NAs...
  if(sum(is.na(data))>0){
    warning("NAs in data => data=na.omit(data) passed to gamlss().")
  }
  
  modelList <- list()
  nms <- c()
  if (is.list(formula)){
    for (i in 1:length(formula)){
      nms <- c(nms,all.names(as.formula(as.character(unname(formula[i])))))
    }
  } else{
    nms <- c(nms,all.names(formula))
  }
  data <- as.data.frame(data)
  
  if(parallel){
    
    no_cores <- length(unique(data$kfold))
    cl <- makeCluster(no_cores)
    registerDoSNOW(cl)
    iterations <- length(unique(data$kfold))
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    modelList <- foreach(fold = unique(data$kfold),.packages = c("gamboostLSS",pckgs),.options.snow = opts) %dopar% {
      
      temp <- gamboostLSS(data = na.omit(data[data$kfold!="Test" & data$kfold!=fold,
                                              which(colnames(data)%in%c(nms))]),
                          formula = formula,
                          families = families,
                          ...)
      
    }
    close(pb)
    stopCluster(cl)
    names(modelList) <- unique(data$kfold)
    
    return(modelList)
    
  } else{
    ### Training Data: k-fold cross-validation/out-of-sample predictions
    
    
    for(fold in unique(data$kfold)){
      
      temp <- gamboostLSS(data = na.omit(data[data$kfold!="Test" & data$kfold!=fold,
                                              which(colnames(data)%in%c(nms))]),
                          formula = formula,
                          families = families,
                          ...)
      
      modelList[[fold]] <- temp
      
    }
    
    return(modelList)

    
  }
  
  
}


#' Convert gamboostLSS objects to MultiQR, or alternatively return predicted parameters of predictive distribution.
#'
#' @param data A \code{data.frame} containing explanatory variables.
#' @param models A Para_gamboostLSS object.
#' @param quantiles Vector of quantiles to be calculated
#'
#' @details Details go here...
#' @return A \code{MultiQR} object derived from gamlss predictive distributions. Alternatively, a matrix condaining the parameters of the predictive gamboostLSS distributions.
#' @export
gamboostLSS_2_MultiQR <- function(data,models,quantiles=seq(0.05,0.95,by=0.05),params=F){
  
  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    data$kfold<-1
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }
  
  data <- as.data.frame(data)
  # Initialise containers for parameters and quantile forecasts
  parameters <- matrix(1,nrow=nrow(data),ncol=4)
  colnames(parameters) <- c("mu", "sigma", "nu", "tau")
  
  distFamily <- c()
  
  ### current implementation requires all data from gamboost model (don't think you can turn off anyway)
  # models[[fold]]$mu$baselearner$`bbs(SWH)`$get_names() <----possible alternatives.....
  # gamb$mu$baselearner$names(gamb$mu$coef())
  
  
  for(fold in unique(data$kfold)){
    
    tempdata <- data[,which(colnames(data)%in%c(colnames(attributes(models[[fold]])$data),"kfold"))]
    
    # NAs not allowed in newdata. Flags required to record possition.
    gooddata <- rowSums(is.na(tempdata))==0
    
    tempPred <- predict(models[[fold]], newdata = tempdata[data$kfold==fold & gooddata,],type="response")
    
    for(i in 1:(length(tempPred))){
      parameters[data$kfold==fold & gooddata,i] <- tempPred[[i]]
    }
    
    distFamily <- unique(c(distFamily,attributes(attributes(models[[fold]])$families)$name))
    
  }
  
  if(params){
    return(parameters)
  }
  
  
  if(length(distFamily)!=1){stop("length(distFamily)!=1 - Only a single parametric distribution family is allowed.")}
  
  myqfun <- attributes(attributes(models[[fold]])$families)$qfun
  
  input <- list()
  if("mu"%in%names(as.list(args(myqfun)))){
    input$mu=parameters[,1]
  }
  if("sigma"%in%names(as.list(args(myqfun)))){
    input$sigma=parameters[,2]
  }
  if("nu"%in%names(as.list(args(myqfun)))){
    input$nu=parameters[,3]
  }
  if("tau"%in%names(as.list(args(myqfun)))){
    input$tau=parameters[,4]
  }
  
  multipleQuantiles <- matrix(NA,nrow=nrow(data),ncol=length(quantiles))
  
  for(i in 1:length(quantiles)){
    input$p <- quantiles[i]
    multipleQuantiles[,i] <- do.call(myqfun,input)
    
  }
  
  colnames(multipleQuantiles) <- paste0("q",100*quantiles)
  multipleQuantiles <- as.data.frame(multipleQuantiles)
  class(multipleQuantiles) <- c("MultiQR",class(multipleQuantiles))
  
  return(multipleQuantiles)
  
  
}


#' Probability integral transform for Para_gamboostLSS objects
#'
#' This function produces a fan plot of a MultiQR object.
#' @param models A Para_gamboostLSS object.
#' @param data Input data corresponding to \code{qrdata}.
#' @param dist_fun cumulative distribution function corresponging to families specified in gamboostLSS model (see example).
#' @param response_name name of response variable in \code{data} object.
#' @details Details go here...
#' @return The probability integral transform of \code{data} through the predictive distribution defined by a list of gamboostLSS objects.
#' @export
gamboostLSS_2_PIT <- function(models,data,dist_fun,response_name,...){
  
  # Arrange kfold cross-validation
  if(is.null(data$kfold)){
    if(length(models)!=1){stop("kfold inconsistent with ppd.")}
    data$kfold<-names(models)
  }else{
    data$kfold[is.na(data$kfold)] <- "Test"
  }
  
  data <- as.data.frame(data)
  
  distFamily <- c()
  
  for(fold in unique(data$kfold)){
    distFamily <- unique(c(distFamily,attributes(attributes(models[[fold]])$families)$name))
  }
  
  if(length(distFamily)!=1){stop("length(distFamily)!=1 - Only a single parametric distribution family is allowed.")}
  
  parameters <- gamboostLSS_2_MultiQR(data=data,models=models,params=T)
  
  
  tempdata <- data[,which(colnames(data)%in%c(colnames(attributes(models[[fold]])$data)))]
  gooddata <- rowSums(is.na(tempdata))==0
  
  input <- list(q=data[[response_name]][gooddata])
  
  if("mu"%in%names(as.list(args(dist_fun)))){
    input$mu=parameters[gooddata,1]
  }
  if("sigma"%in%names(as.list(args(dist_fun)))){
    input$sigma=parameters[gooddata,2]
  }
  if("nu"%in%names(as.list(args(dist_fun)))){
    input$nu=parameters[gooddata,3]
  }
  if("tau"%in%names(as.list(args(dist_fun)))){
    input$tau=parameters[gooddata,4]
  }
  
  X <- rep(NA,nrow(data))
  X[gooddata] <- do.call(dist_fun,input)
  
  
  return(X)
}



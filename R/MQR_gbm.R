#' Multiple Quantile Regression Using GBM
#'
#' This function fits multiple quantile regreesion GBMs with facilities for cross-validation.
#' @param data A \code{data.frame} containing target and explanatory variables. May optionally contain a collumn called "kfold" with numbered/labeled folds and "Test" for test data.
#' @param formaul A \code{formula} object with the response on the left of an ~ operator, and the terms, separated by + operators, on the right
#' @param quantiles The quantiles to fit models for.
#' @param gbm_params List of parameters to be passed to \code{fit.gbm()}.
#' @param CVfolds Control for cross-validation if not supplied in \code{data}.
#' @param perf.plot Plot GBM performance?
#' @param pred_ntree predict using a user-specified tree, if unspecified an out-of-the bag estimate will be used
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
                    pred_ntree = NULL,
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
      
      # Initiate cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
      #set up progress bar
      iterations <- length(unique(data$kfold))
      pb <- txtProgressBar(max = iterations, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      # fit each quantiel model change parameters for CV results
      
      qpred <- foreach(fold = unique(data$kfold),.packages = c("gbm",pckgs),.options.snow = opts) %dopar% {
        
        ### Fit gbm model
        temp_gbm <- do.call(gbm,c(list(formula=formula,data=data[data$kfold!=fold & data$kfold!="Test" & !is.na(data[[formula[[2]]]]),],distribution = list(name="quantile",alpha=q)),gbm_params))
        
        
        ### Save out-of-sample predictions
        if(is.null(pred_ntree)){
          predict.gbm(temp_gbm,
                      newdata = data[data$kfold==fold,],
                      n.trees = gbm.perf(temp_gbm,plot.it = perf.plot))
        } else{
          
          predict.gbm(temp_gbm,
                      newdata = data[data$kfold==fold,],
                      n.trees = pred_ntree)
          
        }
        
        
      }
      
      close(pb)
      stopCluster(cl)
      
      names(qpred) <- unique(data$kfold)
      
      for(fold in unique(data$kfold)){# Loop over CV folds and test data
        
        predqs[[paste0("q",100*q)]][data$kfold==fold] <- qpred[[fold]]
        
      }
    }
    
    
  } else{
    ### Training Data: k-fold cross-validation/out-of-sample predictions
    for(q in quantiles){ # Loop over quantiles
      for(fold in unique(data$kfold)){# Loop over CV folds and test data
        
        
        ### Fit gbm model
        temp_gbm <- do.call(gbm,c(list(formula=formula,data=data[data$kfold!=fold & data$kfold!="Test" & !is.na(data[[formula[[2]]]]),],distribution = list(name="quantile",alpha=q)),gbm_params))
        
        ### Save out-of-sample predictions
        if(is.null(pred_ntree)){
          predqs[[paste0("q",100*q)]][data$kfold==fold] <- predict.gbm(temp_gbm,
                                                                       newdata = data[data$kfold==fold,],
                                                                       n.trees = gbm.perf(temp_gbm,plot.it = perf.plot))
        } else{
          
          predqs[[paste0("q",100*q)]][data$kfold==fold] <- predict.gbm(temp_gbm,
                                                                       newdata = data[data$kfold==fold,],
                                                                       n.trees = pred_ntree)
          
        }
        
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
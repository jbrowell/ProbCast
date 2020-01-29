#' Multiple Quantile Regression Using qgam
#'
#' This function fits multiple quantile generalized additive models (QGAMs) with facilities for cross-validation.
#' @param data A \code{data.frame} containing target and explanatory variables. May optionally contain a collumn called "kfold" with numbered/labeled folds and "Test" for test data.
#' @param formula A \code{formula} object with the response on the left of an ~ operator, and the terms, separated by + operators, on the right
#' @param quantiles The quantiles to fit models for.
#' @param CVfolds Control for cross-validation if not supplied in \code{data}.
#' @param argGam List of parameters to be passed to internal call to \code{mgcv::gam}. See \code{?qgam} for details.
#' @param parallel \code{boolean} parallelize cross-validation process?
#' @param pckgs if parallel is TRUE then  specify packages required for each worker (e.g. c("data.table) if data stored as such)
#' @param cores if parallel is TRUE then number of available cores
#' @param Sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param SortLimits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param save_models_path Path to save models. Model details and file extension pasted onto this string.
#' @details Actual fitting happens in \code{qgam::qgam}.
#' @return Quantile forecasts in a \code{MultiQR} object.
#' @keywords Quantile Regression
#' @importFrom qgam qgam
#' @export
MQR_qreg_qgam <- 
  function(data,
           formula,
           quantiles=c(0.25,0.5,0.75),
           CVfolds=NULL,
           ...,
           argGam = list(),
           w=rep(1,nrow(data)),
           parallel = F,
           cores = NULL,
           pckgs = NULL,
           Sort=T,SortLimits=NULL,
           save_models_path=NULL){
    
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
    
    if(!"BadData" %in% colnames(data)){
      data[,BadData:=F]
    }
    
    
    ### Creae Container for output
    predqs <- data.frame(matrix(NA,ncol = length(quantiles), nrow = nrow(data)))
    colnames(predqs) <- paste0("q",100*quantiles)
    
    ### Model Fitting and Prediction
    if(parallel){
      
      if(is.null(cores)){cores<-detectCores()-1}
      
      ## Do CV in parallel, loop over quantiles
      for(q in quantiles){
        print(paste0("q",q*100))
        
        
        cl <- makeCluster(cores)
        registerDoSNOW(cl)
        iterations <- length(unique(data$kfold))
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        
        predqs_temp <- foreach(fold = unique(data$kfold),
                               .packages = c("qgam",pckgs),
                               .options.snow = opts) %dopar% {
                                 
                                 argGam$weights <- w[data$kfold!=fold & data$kfold!="Test" & data$BadData==F & !is.na(data[[formula[[2]]]])]
                                 
                                 ### Fit model
                                 temp_model <- 
                                   do.call(qgam,c(list(form=formula,
                                                         data=data[data$kfold!=fold & data$kfold!="Test" & data$BadData==F & !is.na(data[[formula[[2]]]]),],
                                                         qu=q,
                                                         argGam=argGam),
                                                    ...))
                                 
                                 if(!is.null(save_models_path)){
                                   save(temp_model,file = paste0(save_models_path,"_q",q,"_fold",fold,".Rda"))
                                 }
                                 
                                 
                                 ### Return out-of-sample predictions
                                 predict(temp_model,
                                         newdata = data[data$kfold==fold,])
                                 
                                 
                                 
                                 
                               }
        close(pb)
        stopCluster(cl)
        
        # Store predictions
        names(predqs_temp) <- unique(data$kfold)
        for(fold in unique(data$kfold)){
          predqs[[paste0("q",100*q)]][data$kfold==fold] <- predqs_temp[[fold]]
        }
        
      }
      
    }else{
      
      ## Loop over quantiles and CV folds
      for(q in quantiles){ # Loop over quantiles
        for(fold in unique(data$kfold)){# Loop over CV folds and test data
          
          argGam$weights <- w[data$kfold!=fold & data$kfold!="Test" & data$BadData==F & !is.na(data[[formula[[2]]]])]
          
          ### Fit model
          temp_model <- do.call(qgam,c(list(form=formula,
                                              data=data[data$kfold!=fold & data$kfold!="Test" & data$BadData==F & !is.na(data[[formula[[2]]]]),],
                                              qu=q,
                                            argGam=argGam),
                                         ...))
          
          ### Save out-of-sample predictions
          predqs[[paste0("q",100*q)]][data$kfold==fold] <- predict(temp_model,
                                                                   newdata = data[data$kfold==fold,])
          
          
          if(!is.null(save_models_path)){
            save(temp_model,file = paste0(save_models_path,"_q",q,"_fold",fold,".Rda"))
          }
          ## Plot Model Coeffs
          # plot(temp_model, which=1)
          # plot(temp_model, which=5)
          # plot(temp_model, which=12)
          
        }
      }
    }
    
    class(predqs) <- c("MultiQR","data.frame")
    
    
    if(Sort){
      predqs <- SortQuantiles(data = predqs,Limits = SortLimits)
    }
    
    return(predqs)
    
  }

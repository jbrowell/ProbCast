###########################################
## MQR using linear quantile regression ####
###########################################
#' Multiple quantile regression, and quantile dressing, using linear quantile regression
#'
#' This function fits multiple conditional linear quantile regression models, optionally to the residuals of
#' a user-specified deterministic forecast with facilities for cross-validation.
#' 
#' @author Jethro Browell, \email{jethro.browell@@glasgow.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory
#' variables. Optionally, supplying a deterministic forecast with \code{offset} will return
#' the this forecast dressed with multiple predictive quantiles.
#' May optionally contain a column called "kfold" with numbered/labelled folds and "Test" for test data. 
#' @param formala A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right passed to \code{rq()}.
#' @param quantiles The quantiles to fit models for.
#' @param offset The column name in \code{data} of an optional deterministic forecast to be dressed with multiple predictive quantiles.
#' @param cv_folds Control for cross-validation with various options, either:
#' \itemize{
#'  \item the column name of the fold index supplied in data. Observations and inputs 
#'   in the index labelled "Test" will serve as test data and held out in model training.
#'  \item an integer giving the number of cross validation folds to generate. Folds are constructed as block chunks. 
#'  Default behaviour is 5 folds.
#'  \item NULL indicates that no cross validation should be performed and the returned model is trained on all \code{data}.
#' }
#' @param exclude_train A column name in \code{data} indicating if a row should be excluded from model
#' training, i.e. if it contains bad data (will be coerced to \code{logical}). Alternatively,
#' an \code{integer} or \code{logical} vector with length equal to the number of rows in \code{data} indicating
#' the same. Rows labelled \code{TRUE} are excluded from model training.
#' @param sort \code{boolean} Sort quantiles using \code{SortQuantiles()}?
#' @param sort_limits \code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.
#' @param ... Additional arguments passed to \code{rq}.
#' @details The returned predictive quantiles are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' @return Returns a \code{list} containing predictive quantiles (in a \code{MultiQR} object) and \code{rq} models.
#' @keywords Quantile Regression
#' @export
qreg_mrq <- function(data,
                     formula,
                     quantiles=c(0.25,0.5,0.75),
                     offset = NULL,
                     cv_folds=NULL,
                     exclude_train=NULL,
                     sort=T,sort_limits=NULL,
                     ...){
  
  
  # Set-up cv folds & do checks
  cv_labs <- cv_control(data = data,cv_folds = cv_folds)
  
  # Exclude points from training? & do checks
  exclude_idx <- exclude_fun(data = data,exclude_train = exclude_train)
  
  # Initialise output container
  FINAL_OUTPUT <- list(mqr_pred=data.table(matrix(as.numeric(NA),
                                                  ncol = length(quantiles),
                                                  nrow = nrow(data))),
                       call=match.call(),
                       kfold_index = cv_labs$idx,
                       model_names = cv_labs$fold_loop,
                       exclude_index = exclude_idx,
                       models = list(),
                       sorted = list(sort=sort,
                                     sort_limits=sort_limits),
                       default_model=if("Test"%in%cv_labs$fold_loop){"Test"}else{cv_labs$fold_loop[1]})
  colnames(FINAL_OUTPUT$mqr_pred) <- paste0("q",100*quantiles)
  class(FINAL_OUTPUT) <- c("qreg_mrq",class(FINAL_OUTPUT))
  
  # Adjust formula if offset being used (rq does not handle offsets automatically in the same way as lm)
  if(!is.null(offset)){
    formula = as.formula(paste0(c(terms(formula)[[2]],"-",offset," ~ ",terms(formula)[[3]]),collapse = ""))
  }
  
  for(fold in FINAL_OUTPUT$model_names){
    
    print(paste0("Multiple rq, kfold=",fold))
    
    for(i in 1:length(quantiles)){
      
      idx <- FINAL_OUTPUT$kfold_index!=fold & FINAL_OUTPUT$kfold_index!="Test" & FINAL_OUTPUT$exclude_index==0
      
      rq_model <- rq(formula,
                     tau = quantiles[i],
                     data = data[idx,],
                     ...)
      
      ## Throw away unneeded entries in "rq" as take up a lot of space!!!
      rq_model <- rq_model[c("coefficients","terms","xlevels","contrasts")]
      class(rq_model) <- "rq"
      
      FINAL_OUTPUT$models$rqs[[fold]][[paste0("q",100*quantiles[i])]] <- rq_model
      
      if(is.null(cv_folds)){
        FINAL_OUTPUT$mqr_pred[,i] <- predict.rq(FINAL_OUTPUT$models$rqs[[fold]][[paste0("q",100*quantiles[i])]],data)
      }else{
        FINAL_OUTPUT$mqr_pred[FINAL_OUTPUT$kfold_index==fold,i] <- predict.rq(FINAL_OUTPUT$models$rqs[[fold]][[paste0("q",100*quantiles[i])]],data[FINAL_OUTPUT$kfold_index==fold,])
      }
    }
  }
  
  
  ## Order columns
  FINAL_OUTPUT$mqr_pred <- FINAL_OUTPUT$mqr_pred[,order(as.numeric(gsub("q","",colnames(FINAL_OUTPUT$mqr_pred)))),with=F]
  
  ## Add back offset
  if(!is.null(offset)){
    FINAL_OUTPUT$mqr_pred <- FINAL_OUTPUT$mqr_pred + data[[offset]]
  }
  
  ## Assign MultiQR class
  if(class(FINAL_OUTPUT$mqr_pred)[1]!="MultiQR"){
    class(FINAL_OUTPUT$mqr_pred) <- c("MultiQR",class(FINAL_OUTPUT$mqr_pred))
  }
  
  ## Sort quantiles
  if(sort){
    FINAL_OUTPUT$mqr_pred <- SortQuantiles(data = FINAL_OUTPUT$mqr_pred,
                                           Limits = sort_limits)
  }
  
  return(FINAL_OUTPUT)
}

#' cross validation set-up for model fitting functions
#'
#' This generates kfold indexes and vectors of reference values for the model fitting loop  
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory
#' variables. May optionally contain a column with labeled folds and "Test" for test data. See \code{cv_folds}.
#' @param cv_folds Control for cross-validation with various options, either:
#' \itemize{
#'  \item the column name of the fold index supplied in data. Using this option, if there is a 
#'   fold called "Test", this will serve as test data and held out in the model training.
#'  \item an integer giving the number of cv evaluations to perform. Folds are constructed as block chunks. 
#'  Default behaviour is 5 folds.
#'  \item NULL indicates that no cross validation should be performed and the returned model is trained on all \code{data}.
#' }
#' @details ....
#' @return A named list containing the kfold indexes and the kfold loop reference values
cv_control <- function(data,
                       cv_folds = "kfold"){
  
  
  ## - cv from data
  if(is.character(cv_folds)){
    if(sum(cv_folds %in% colnames(data))==0){
      
      stop(paste0("cannot find column ",cv_folds," in data.") )
           
    }
    
    fold_idx <- data[[cv_folds]]
    
    if(is.numeric(fold_idx)){
      
      fold_idx <- paste0("fold",fold_idx)
      
    }
    fold_loop <- unique(fold_idx)
    
    if(length(fold_loop)<=1){
      
      stop("number of unique kfolds in data should be > 1. Alternatively, set cv_folds to NULL for no cv operation")
    }
    
    if(sum(is.na(fold_idx))>0){
      
      stop(paste0("NAs detected in data$",cv_folds))
      
    }
    
    
    
  }
  
  ## -  no cv
  if(is.null(cv_folds)){
    
    
    fold_idx <- rep(1,nrow(data))
    fold_loop <- "all_data"
    
    
  }
  
  ## - auto cv
  if(is.numeric(cv_folds)){
    
    if(cv_folds<=1){
      
      stop("cv_folds should be > 1. Alternatively, set cv_folds to NULL for no cv operation") 
    }
    
    cv_folds <- as.integer(cv_folds)
    fold_idx <- paste0("fold",sort(rep(1:cv_folds,length.out=nrow(data))))
    fold_loop <- paste0("fold",1:cv_folds)
    
    
  }
  
  
  return(list(idx=fold_idx,fold_loop=fold_loop))
  
}
  
  
  
  
  
  

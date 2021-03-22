
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
#'  \item vector of length==nrow(data) containing character or numeric fold labels
#' }
#' @details ....
#' @return A named list containing the kfold indexes and the kfold loop reference values
cv_control <- function(data,
                       cv_folds = "kfold"){
  
  ## -  no cv
  if(is.null(cv_folds)){
    
    
        fold_idx <- rep(1,nrow(data))
        fold_loop <- "all_data"
    
  ## - cv folds character(data) or integer(auto_cv)  
  } else if(length(cv_folds)==1){
    
        ## - cv from data
        if(is.character(cv_folds)){
          if(sum(cv_folds %in% colnames(data))==0){
            
            stop(paste0("cannot find column ",cv_folds," in data.") )
          }
          
          fold_idx <- data[[cv_folds]]
          if(is.numeric(fold_idx)){
            fold_idx <- paste0("fold",fold_idx)
          }
        
          # else auto-cv  
        } else  if(is.numeric(cv_folds)){
          
          if(cv_folds<=1){
            
            stop("cv_folds should be > 1. Alternatively, set cv_folds to NULL for no cv operation") 
          }
          
          cv_folds <- as.integer(cv_folds)
          fold_idx <- paste0("fold",sort(rep(1:cv_folds,length.out=nrow(data))))
          
          
        } else {stop("check cv_folds argument")}
        
        fold_loop <- unique(fold_idx)
    
    
    # cv folds user-supplied direct vector
  } else if(length(cv_folds)==nrow(data)){
    
    
      if(is.character(cv_folds)){
        
        fold_idx <- cv_folds
        
      } else  if(is.numeric(cv_folds)){
        
        ## else change to character
        fold_idx <- paste0("fold",cv_folds)
        
      } else {stop("check cv_folds argument")}
      
      
      fold_loop <- unique(fold_idx)
    
    
  } else{stop("check cv_folds argument")}
  
  
  # checks
  if(!("all_data"%in%fold_loop) & length(fold_loop)<=1){
    
    stop("number of unique kfolds in data should be > 1. Alternatively, set cv_folds to NULL for no cv operation")
  }
  
  if(sum(is.na(fold_idx))>0){
    
    stop("NAs detected in cv_folds")
    
  }
  
  
  # for kfold idx is always a character vector, except when cv_folds=NULL
  # fold_loop is always a character
  return(list(idx=fold_idx,fold_loop=fold_loop))
  
}







#'  set up exclude index for modelling fitting
#'
#' This generates a vector of length nrow(data) to flag removal from model **training**
#' 
#' @author Ciaran Gilbert, \email{ciaran.gilbert@@strath.ac.uk}
#' @param data A \code{data.frame} containing target and explanatory
#' variables. May optionally contain a column of exclude flags named by \code{exclude_train}.
#' @param exclude_train control for exclusion of rows in data for the model training only, with various options, either:
#' \itemize{
#'  \item the column name of the binary/boolean exclude flag supplied in data.
#'  \item a vector of binary/boolean exclusion flags of length nrow(data)
#'  \item NULL indicates no exclusion
#' }
#' @details ....
#' @return A vector of binary flags for removing training data points
exclude_fun <- function(data,exclude_train = NULL){
  
  
  # exclude from data
  if(is.character(exclude_train)){
    
    if(sum(exclude_train %in% colnames(data))==0){
      stop(paste0("cannot find column ",exclude," in data."))
    }
    
    exclude_idx <- as.numeric(data[[exclude_train]])
    

  # exclude from vector  
  } else if(!is.null(exclude_train)){
    
    if(nrow(data)!=length(exclude_train)){
      stop("nrow(data)!=length(exclude_train)")
    }
    
    exclude_idx <- as.numeric(exclude_train)
    
  # no exclusion  
  } else{
    
    exclude_idx <- rep(0,nrow(data))
    
  }
  
  # if !(indx%in%c(0,1))
  if(sum(Negate('%in%')(unique(exclude_idx),c(0,1)))!=0){
    
    stop(paste0("exclude_train should be a binary or boolean vector"))
    
  }
  
  
  if(sum(is.na(exclude_idx))>0){
    
    stop(paste0("NAs detected in exlude_train"))
    
  }
  
  return(exclude_idx) 
  
}
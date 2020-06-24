#' Extreme value distirbutions for distribution tails
#'
#' This function fits semi-parametric models for extreme value distributions
#' using \code{evgam} with facilities for cross-validation.
#' 
#' @author Jethro Browell, \email{jethro.browell@@strath.ac.uk}
#' @param data The modelling table as a \code{data.frame} or \code{data.table}
#' @param mqr_data A \code{MultiQR} object corresponding to \code{data}
#' @param tail_starts The upper and lower probability levels beyond which the
#'  tail distribution is used. E.g. \code{c(5,95)} to use tail distribution for probabilities
#'  below quantile 0.05 and above 0.95. Defaults to min and max from \code{mqr_data} 
#' @param formala A \code{formula} object with the response on the left
#' of an ~ operator, and the terms, separated by + operators, on the right
#' @param CVfolds Control for cross-validation if not supplied in \code{data}
#' @param BadData_col Name of a boolean column in \code{data} indicating bad data
#' which should be excluded from parameter estimation
#' 
#' @details The returned predictive quantiles are those produced out-of-sample for each
#' cross-validation fold (using models trained on the remaining folds but not "Test" data).
#' Predictive quantiles corresponding to "Test" data are produced using models trained on all
#' non-test data.
#' @return Returns \code{data} with additional columns containing the predicted
#' parameters of the specified tail distirbution
#' @keywords Extreme Value Distribution; Tails
#' @import evgam
#' @import data.table
#' @export
tails_ev <- function(data,
                     mqr_data,
                     tail_starts=range(as.numeric(gsub("q","",names(mqr_data))),na.rm = T),
                     formula,
                     CVfolds=NULL,
                     BadData_col=NULL){
  
  ### Add BadData column if doesn't exist
  if(!exists(data$BadData) & is.null(BadData_col)){
    ## No bad data
    data[,BadData:=F]  
  }else{
    ## BadData indicator from a different column
    data[,BadData:=get(BadData_col)]  
  }
  
  ### Check rows of data and mqr_data
  if(nrow(mqr_data)!=nrow(data)){stop("nrow(mqr_data) must equal nrow(data).")}
  
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
  
  ### Convert input data to data.table
  data <- as.data.table(data)
  
  ### Get target variable name from formula
  target <- "???"
  
  
  ### Get "tail residuals"
  data$tail_l_resid <- -(data[,get(target)] - mqr_data[[paste0("q",tail_starts[1])]])
  data$tail_r_resid <- data[,get(target)] - mqr_data[[paste0("q",tail_starts[2])]]
  
  
  ### Fit model in CV
  stop("Function not finished!!!")
  for(fold in unique(data$kfold)){# Loop over CV folds and test data
    
    fit_l <- evgam(list(tail_l_resid~1,
                        ~1),
                   data = data[tail_l_resid>0 & BadData==F & !(kfold%in%c(fold,"Test")),],
                   family = "gpd")
    
    fit_r <- evgam(list(tail_r_resid~1,
                        ~1), data = data[tail_r_resid>0 & BadData==F & !(kfold%in%c(fold,"Test")),],
                   family = "gpd")
    
    # summary(fit_l)
    # plot(fit_l)    
    
    
    ## Estimate Parameters
    gpd_params_l <- predict(fit_l,newdata = NodeData[[N]])
    gpd_params_r <- predict(fit_r,newdata = NodeData[[N]])
  }
  
  return(data)
  
  
}

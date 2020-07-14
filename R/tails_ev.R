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
#' of an ~ operator, and the terms, separated by + operators, on the right. This formula is
#' used for both upper and lower tails unless \code{formula_r} is specified in which case
#' \code{formula} is used for the lower tail only.
#' @param formula_r as \code{formula} but for the upper tail only.
#' @param CVfolds Control for cross-validation if not supplied in \code{data}.
#' @param BadData_col Name of a boolean column in \code{data} indicating bad data
#' which should be excluded from parameter estimation.
#' @param evgam_family specifies distribution, see \code{?evgam}.
#' @param print_summary If \code{TRUE}, summary is printed for each evgam fit.
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
                     formula_r=formula,
                     CVfolds=NULL,
                     BadData_col=NULL,
                     evgam_family = "gpd",
                     print_summary=F){
  
  ## Input Checks
  if(evgam_family!="gpd"){warning("Only tested for evgam_family = \"gpd\"...")}
  if(length(tail_starts)!=2){stop("tail_starts should be length 2 and in paste0(\"q\",names(mqr_data)).")}
  
  ### Add BadData column if doesn't exist
  if(is.null(data$BadData) & is.null(BadData_col)){
    ## No bad data
    data[,BadData:=F]  
  }else if(is.null(data$BadData)){
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
  
  ### Get target variable name from formula and prep formulas for evgam
  target <- paste0(formula[[1]][[2]])
  formula[[1]] <- reformulate(deparse(formula[[1]]),response="tail_l_resid")
  formula_r[[1]] <- reformulate(deparse(formula_r[[1]]),response="tail_r_resid")
  
  ### Get "tail residuals"
  data$tail_l_resid <- -(data[,get(target)] - mqr_data[[paste0("q",tail_starts[1])]])
  data$tail_r_resid <- data[,get(target)] - mqr_data[[paste0("q",tail_starts[2])]]
  
  
  ### Fit model using CV
  for(fold in unique(data$kfold)){# Loop over CV folds and test data
    
    fit_l <- evgam(formula,
                   data = data[tail_l_resid>0 & BadData==F & !(kfold%in%c(fold,"Test")),],
                   family = evgam_family)
    
    fit_r <- evgam(formula_r,
                   data = data[tail_r_resid>0 & BadData==F & !(kfold%in%c(fold,"Test")),],
                   family = evgam_family)
    
    if(print_summary){
      print(summary(fit_l))
      print(summary(fit_r))
      
      plot(fit_l) 
      plot(fit_r) 
    }
    
    
    ## Estimate Parameters
    params_l <- predict(fit_l,newdata = data[kfold==fold],type="response")
    params_r <- predict(fit_r,newdata = data[kfold==fold],type = "response")
    
    data[kfold==fold,c(paste0(evgam_family,"_",names(params_l),"_l"),
                       paste0(evgam_family,"_",names(params_r),"_r")):=cbind(params_l,params_r)]
    
  }
  
  return(data)
  
  
}

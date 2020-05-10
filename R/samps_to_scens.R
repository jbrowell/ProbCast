#' Generate multivariate forecasts
#'
#' This function produces a list of multivariate scenario forecasts in the marginal domain from the spatial/tempral/spatiotemporal gaussian covariance matrices and marginal distributions
#' @param copulatype Either "spatial" of "temporal", note that spatio-temporal can be generated via "temporal" setting
#' @param no_samps Number of scenarios required
#' @param marginals a named list of the margins of the copula - e.g. if class is MultiQR --> list(<<name>> = <<MultiQR object>>). Multiple margins are possible for multiple locations (see examples) although they must be the same class (MQR or distribution parameters). If parametric class supply a list of the distribution parameters here and the corresponding quantile function in \code{control} (see below). The ordering of this list is important for multiple locations --- it should be ordered according to the row/columns in each member of \code{sigma_kf}
#' @param sigma_kf a named list of the covariance matrices with elements corresponding to each fold.
#' @param mean_kf a named list of the mean vectors with elements corresponding to each fold
#' @param control a named list of with nested control parameters (named according to \code{marginals}). Each named list should contain \code{kfold}, \code{issue_ind}, and \code{horiz_ind} which are the kfold, issue time, and lead time vectors corresponding to the margins of the copula. If margins are MultiQR class also define \code{PIT_method} and list \code{CDFtails}, which are passed to the PIT function. If the margins are distribution parameter predictions then define \code{q_fun}, which transforms the columns of \code{marginals} through the quantile function --- see example for more details. 
#' @param mcmapply_cores defaults to 1. Warning, only change if not using windows OS --- see the \code{parallel::mcmapply} help page for more info. Speed improvements possible when generating sptio-temporal scenarios, set to the number of locations if possible.
#' @param mvnfast_cores defaults to 1. See \code{mvnfast::rmvn}
#' @param ... other parameters to be passed to mvtnorm::rmvnorm
#' @note For spatio-temporal scenarios, each site must have the same number of inputs to the governing covariance matrix.
#' @note For multiple locations the ordering of the lists of the margins & control, and the structure of the covariance matrices is very important; if the columns/rows in each covariance matrix are ordered loc1_h1, loc1_h2,..., loc2_h1, loc2_h_2,..., loc_3_h1, loc_3_h2,... i.e. location_leadtime --- then the list of the marginals should be in the same order loc1, loc2, loc3,....
#' @note Ensure kfold ids in the control list do not change within any issue time --- i.e. make sure the issue times are unique to each fold. 
#' @details Details go here...
#' @return A list of \code{data.table} objects containing multivariate scenario forecasts
#' @examples
#' \dontrun{
#' # for parametric type marginals with a Generalized Beta type 2 family
#' scens <- samps_to_scens(copulatype = "temporal",no_samps = 100,marginals = list(loc_1 = param_margins),sigma_kf = cvm,mean_kf = mean_vec,
#'                         control=list(loc_1 = list(kfold = loc_1data$kfold, issue_ind = loc_1data$issue_time, horiz_ind = loc_1data$lead_time,
#'                                                   q_fun = gamlss.dist::qGB2)))
#' }
#' \dontrun{
#' # for MQR type marginals
#' scens <- samps_to_scens(copulatype = "temporal",no_samps = 100,marginals = list(loc_1 = mqr_gbm_1),sigma_kf = cvm,mean_kf = mean_vec,
#'                         control=list(loc_1 = list(kfold = loc_1data$kfold, issue_ind = loc_1data$issue_time, horiz_ind = loc_1data$lead_time,
#'                                                   PIT_method = "linear",CDFtails= list(method = "interpolate", L=0,U=1))))
#' }
#' \dontrun{
#' # for spatio-temporal scenarios with MQR type marginals
#' scens <- samps_to_scens(copulatype = "temporal", no_samps = 100,marginals = list(loc_1 = mqr_gbm,loc_2 = mqr_gbm_2),sigma_kf = cvm_2,mean_kf = mean_vec_2,
#'                         control=list(loc_1 = list(kfold = loc_1data$kfold, issue_ind = loc_1data$issue_time, horiz_ind = loc_1data$lead_time,
#'                                                   PIT_method = "linear",CDFtails= list(method = "interpolate", L=0,U=1)),
#'                                      loc_2 = list(kfold = loc_2data$kfold, issue_ind = loc_2data$issue_time, horiz_ind = loc_2data$lead_time,
#'                                                   PIT_method = "linear", CDFtails = list(method = "interpolate", L=0, U=1))))
#' }
#' @importFrom mvnfast rmvn
#' @importFrom parallel mcmapply
#' @import data.table
#' @export
samps_to_scens <- function(copulatype,no_samps,marginals,sigma_kf,mean_kf,control,mcmapply_cores = 1L, mvnfast_cores = 1L,...){
  
  # no kfold capability?
  # improve ordering of lists cvm matrix...
  
  # mean_list <- list()
  # for (i in levels(unique(u_obsind$kfold))){
  #   mean_list[[i]] <- rep(0, 24)
  # }
  # 
  # 
  # copulatype <- "temporal"
  # no_samps <- 200
  # marginals <- list(loc_1 = test1$gbm_mqr)
  # sigma_kf <- cvm_gbm
  # mean_kf <- mean_list
  # control <- list(loc_1 = list(kfold = u_obsind$kfold,issue_ind=u_obsind$i_time,horiz_ind=u_obsind$lead_time,
  #                           PIT_method="spline",
  #                           CDFtails = list(method="interpolate",L=0,U=1,ntailpoints=100)))
  # mcmapply_cores <- 1L
  # mvnfast_cores <- 1L
  
  if(class(marginals)[1]!="list"){
    marginals <- list(loc_1 = marginals)
    warning("1 location detected --- margin coerced to list")
    if(length(control)>1){
      control <- list(loc_1 = control)
    }
    
  }
  
  if(length(marginals)!=length(control)){
    stop("control dimensions must equal marginals")
  }
  
  if(!identical(names(marginals),names(control))){
    stop("control must be named and in the same order as marginals")
  }
  
  if(!identical(names(sigma_kf),names(mean_kf))){
    stop("mean_kf order must equal sigma_kf")
  }
  
  if(mcmapply_cores!=1){
    warning("Only change mcmapply_cores if not using Windows OS")
  }
  
  
  if(copulatype=="spatial"){
    
    
    # This function is for extracting spatial scenario samples
    extr_kf_spatsamp <- function(kf_samp_df,uni_kfold){
      
      print(paste0("taking samples for --- ",uni_kfold))
      
      ## take cholesky decomp. so we only have to do it once for all issuetimes in fold
      chol_mat <- chol(sigma_kf[[uni_kfold]])
      
      arg_list <- list(n=no_samps,sigma=chol_mat,mu=mean_kf[[uni_kfold]], isChol = TRUE, ncores = mvnfast_cores)
      
      # sample from multivariate gaussian, gives results in a list of matrices
      kf_samps <- replicate(n=nrow(kf_samp_df),expr=do.call(eval(parse(text="rmvn")),args=arg_list),simplify = F)
      
      # transform sample rows ---> samples in cols and time_ind in rows
      kf_samps <- lapply(kf_samps,t)
      
      # convert to uniform domain
      kf_samps <- lapply(kf_samps,pnorm)
      
      # add ID row column for split list later (will help margins length>1) --- poss imp, impose naming convention on cvms?
      kf_samps <- lapply(kf_samps,function(x){cbind(x,sort(rep(1:length(marginals),nrow(x)/length(marginals))))})
      
      #bind the rows
      kf_samps <- data.table(docall("rbind",kf_samps))
      
      # split the data.table up into a list of data.tables by the rowID column for different locations
      kf_samps <- split(kf_samps,by=tail(colnames(kf_samps),1),keep.by = FALSE)
      
      # bind with kf_samp_df time indices
      kf_samps <- lapply(kf_samps,function(x){cbind(kf_samp_df,x)})
      
      # name list
      names(kf_samps) <- names(marginals)
      
      return(kf_samps)
      
    }
    
    # find the unique combinations of issue_time and horizon at per fold across all the locations
    find_nsamp <- list()
    for(i in names(sigma_kf)){
      find_nsamp[[i]] <- unique(do.call(rbind,unname(lapply(control,function(x){data.table(issue_ind=x$issue_ind[x$kfold==i],horiz_ind=x$horiz_ind[x$kfold==i])}))))
      find_nsamp[[i]] <- find_nsamp[[i]][order(find_nsamp[[i]]$issue_ind, find_nsamp[[i]]$horiz_ind),]
    }
    
    # extract samples calling etr_kf_spatsamp
    samps <- mapply(extr_kf_spatsamp,kf_samp_df=find_nsamp,uni_kfold = as.list(names(find_nsamp)),SIMPLIFY = F)
    rm(find_nsamp)
    
    
    
  } else{ if (copulatype=="temporal"){
    
    # This function is for extracting temporal/spatio-temporal scenario samples
    extr_kf_temposamp <- function(issuetimes,uni_kfold){
      
      print(paste0("taking samples for --- ",uni_kfold))
      
      ## take cholesky decomp. so we only have to do it once for all issuetimes in fold
      chol_mat <- chol(sigma_kf[[uni_kfold]])
      
      arg_list <- list(n=no_samps,sigma=chol_mat,mu=mean_kf[[uni_kfold]], isChol = TRUE, ncores = mvnfast_cores)
      
      # sample from multivariate gaussian, gives results in a list of matrices
      kf_samps <- replicate(n=length(issuetimes$issue_ind),expr=do.call(eval(parse(text="rmvn")),args=arg_list),simplify = F)
      
      # transform sample rows ---> samples in cols and horizon in rows
      kf_samps <- lapply(kf_samps,t)
      
      # convert to uniform domain
      kf_samps <- lapply(kf_samps,pnorm)
      
      # add ID row column for split list later (will help margins length>1) --- poss imp, impose naming convention on cvms?(and kpnames from mvnfast)
      kf_samps <- lapply(kf_samps,function(x){cbind(x,sort(rep(1:length(marginals),nrow(x)/length(marginals))))})
      
      # bind the rows
      kf_samps <- data.table(docall("rbind",kf_samps))
      
      # add issueTime ID to each data.table
      issue_ind <- sort(rep(issuetimes$issue_ind,nrow(kf_samps)/length(issuetimes$issue_ind)))
      
      # add horiz_ind to vector --- need to change this**** impose naming convention on cvms?(and kpnames from mvnfast)
      horiz_ind <- rep(sort(unique(control[[1]]$horiz_ind)),nrow(kf_samps)/length(sort(unique(control[[1]]$horiz_ind))))
      kf_samps <- cbind(issue_ind,horiz_ind,kf_samps)
      
      # split the data.table up into a list of data.tables by the rowID column for different locations
      kf_samps <- split(kf_samps,by=tail(colnames(kf_samps),1),keep.by = FALSE)
      
      # name list
      names(kf_samps) <- names(marginals)
      
      
      return(kf_samps)
      
    }
    
    
    
    # find number of unique issue_times per fold across all the locations (use data.frame to avoid losing posixct class if present)
    find_issue <- list()
    for(i in names(sigma_kf)){
      find_issue[[i]] <- unique(do.call(rbind,unname(lapply(control,function(x){data.table(issue_ind=unique(x$issue_ind[x$kfold==i]))}))))
      find_issue[[i]] <- find_issue[[i]][order(find_issue[[i]]$issue_ind),]
    }
    
    
    # extract samples calling etr_kf_temposamp for each fold
    samps <- mapply(extr_kf_temposamp,issuetimes=find_issue,uni_kfold = as.list(names(find_issue)),SIMPLIFY = F)
    rm(find_issue)
    
    
  }else{stop("copula type mis-specified")}}
  
  
  
  # merge samples with control data to filter to the required samples for passing to the PIT
  # output from above is samps$<<fold>>$<<loc>> need to rearrage into samps$<<loc>>
  samps <- lapply(names(marginals),function(i){rbindlist(lapply(samps,function(x){x[[i]]}))})
  names(samps) <- names(marginals)
  
  # set order of df
  samps <- lapply(samps,function(x){setorder(x,issue_ind,horiz_ind)})
  # merge control cols and the samples to give the final samps for transormation through PIT
  cont_ids <- lapply(control,function(x){data.table(issue_ind=x$issue_ind,horiz_ind=x$horiz_ind,sort_ind=1:length(x$issue_ind))})
  samps <- mapply(merge.data.table,x = cont_ids,y = samps,MoreArgs = list(all.x=T),SIMPLIFY = F)
  rm(cont_ids)
  ## preserve order of merged scenario table with input control table
  samps <- lapply(samps,function(x){setorder(x,sort_ind)})
  
  
  print(paste0("Transforming samples into original domain"))
  ### transform Unifrom Variable into original domain
  ### add S3 support for PPD...
  if (class(marginals[[1]])[1]%in%c("MultiQR")){
    
    method_list <- lapply(control,function(x){x$PIT_method})
    CDFtail_list <- lapply(control,function(x){x$CDFtails})
    
    marg_names <- lapply(marginals,colnames)
    samp_names <- lapply(samps,function(x){colnames(x)[-c(1:3)]})
    
    ## reduce memory usage in parallel pass_invcdf by joining outside and try gc()
    samps <- mapply(cbind, samps, marginals, SIMPLIFY = F)
    rm(marginals)
    invisible(gc())
    
    ### faster inverse pit using data.table (x7 in example.R)
    pass_invcdf <- function(dt, s_nms, q_nms,...){
      
      ## using as.list here makes data.table return in 'wide format automatically & as numeric(.SD) prevents input samps from being a list..
      dt <- dt[,as.list(contCDF(quantiles = t(mget(q_nms)),inverse = TRUE,...)(as.numeric(.SD))),by=.(sort_ind),.SDcols=s_nms]
      setorder(dt,sort_ind)
      dt[,sort_ind := NULL]
      return(dt)
      
    }
    
    samps <- mcmapply(function(...){pass_invcdf(...)},dt = samps, s_nms = samp_names, q_nms = marg_names, method = method_list,
                           tails = CDFtail_list, SIMPLIFY = F,mc.cores = mcmapply_cores)
    
    
  } else {
    
    ### make sure before as.list in do.call the object is a data.frame
    return_ppdsamps <- function(samps,margin,quant_f){
      ppd_samps <- as.data.table(lapply(samps,function(x){do.call(quant_f,as.list(data.table(cbind(p=x,margin))))}))
      return(ppd_samps)
    }
    
    q_list <- lapply(control,function(x){x$q_fun})
    
    # remove issuetime, horizon, and sorting column for passing through PIT-
    samps <- lapply(samps,function(x){x[,-c(1:3)]})
    samps <- mcmapply(return_ppdsamps,samps = samps, margin = marginals,quant_f = q_list,SIMPLIFY = F,mc.cores = mcmapply_cores)
    
    
  }
  
  samps <- lapply(samps,function(x){`colnames<-`(x,paste0("scen_",1:ncol(x)))})
  
  
  return(samps)
  
  
}

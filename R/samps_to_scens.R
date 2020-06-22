#' Generate multivariate forecasts
#'
#' This function produces a list of multivariate scenario forecasts in the marginal domain from the spatial/tempral/spatiotemporal gaussian covariance matrices and marginal distributions
#' @param copulatype Either "spatial" of "temporal", note that spatio-temporal can be generated via "temporal" setting
#' @param no_samps Number of scenarios required
#' @param marginals a named list of the margins of the copula - e.g. if class is MultiQR --> list(<<name>> = <<MultiQR object>>). Multiple margins are possible for multiple locations (see examples) although they must be the same class (MQR or distribution parameters). If parametric class supply a list of the distribution parameters here and the corresponding quantile function in \code{control} (see below). The ordering of this list is important for multiple locations --- see note below.
#' @param sigma_kf a named list of the covariance matrices with elements corresponding to each fold.
#' @param mean_kf a named list of the mean vectors with elements corresponding to each fold
#' @param control a named list of with nested control parameters (named according to \code{marginals}). Each named list should contain \code{kfold}, \code{issue_ind}, and \code{horiz_ind} which are the kfold, issue time, and lead time vectors corresponding to the margins of the copula. If margins are MultiQR class also define \code{PIT_method} and list \code{CDFtails}, which are passed to the PIT function. If the margins are distribution parameter predictions then define \code{q_fun}, which transforms the columns of \code{marginals} through the quantile function --- see example for more details. 
#' @param mcmapply_cores defaults to 1. Warning, only change if not using windows OS --- see the \code{parallel::mcmapply} for more info. Speed improvements possible when generating sptio-temporal scenarios, set to the number of locations if possible.
#' @param mvnfast_cores defaults to 1. See \code{mvnfast::rmvn}
#' @param chunk_dir a character string containing a directory for storing temporary chunked datasets per fold. Useful for very high dimensional distributions or many samples. Defaults to \code{NULL}, i.e. no chunk. Only valid for \code{Multi.QR} type marginals for now...
#' @param ... extra arguments to \code{contCDF} to be applied to all marginals, e.g to deal with ties.
#' @note For spatio-temporal scenarios, each site must have the same number of lead-times in the covariance matrix.
#' @note For multiple locations the ordering of the lists of the margins and the structure of the covariance matrices is very important; if the columns/rows in each covariance matrix are ordered loc1_h1, loc1_h2,..., loc2_h1, loc2_h_2,..., loc_3_h1, loc_3_h2,... i.e. location_leadtime --- then the list of the marginals should be in the same order loc1, loc2, loc3,....
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
#' @importFrom fst write_fst read_fst
#' @import data.table
#' @export
samps_to_scens <- function(copulatype,no_samps,marginals,sigma_kf,mean_kf,control,mcmapply_cores = 1L, mvnfast_cores = 1L, chunk_dir = NULL,...){
  
  # no kfold capability?
  # improve ordering of lists cvm matrix...
  # imposing a class on cvms will help a lot
  # at the moment the location/leadtime is inferred from data assuming a v. specific ordering of the matrix/marginals
  
  # mean_list <- list()
  # for (i in levels(unique(u_obsind$kfold))){
  #   mean_list[[i]] <- rep(0, 24)
  # }
  # 
  # 
  # copulatype <- "temporal"
  # no_samps <- 2000
  # marginals <- list(loc_1 = test1$gbm_mqr)
  # sigma_kf <- cvm_gbm
  # mean_kf <- mean_list
  # control <- list(loc_1 = list(kfold = u_obsind$kfold,issue_ind=u_obsind$i_time,horiz_ind=u_obsind$lead_time,
  #                           PIT_method="spline",
  #                           CDFtails = list(method="interpolate",L=0,U=1,ntailpoints=100)))
  # mcmapply_cores <- 1L
  # mvnfast_cores <- 1L
  # chunk_dir <- "/Users/Ciaran/"
  # rm(control,loc_list,marginals,method_list,mean_kf,samps,samps_temp,sigma_kf,chunk_dir,copulatype,mcmapply_cores,mvnfast_cores,i,CDFtail_list)
  
  if(class(marginals)[1]!="list"){
    marginals <- list(loc_1 = marginals)
    warning("1 location detected --- margin coerced to list")
    if(length(control)>1){
      control <- list(loc_1 = control)
    }
    
  }
  
  if(length(unique(sapply(marginals,function(x){class(x)[1]})))!=1){
    stop("marginals must all be the same class")
  }
  
  if(class(marginals[[1]])[1]!="MultiQR" & !is.null(chunk_dir)){
    warning("chunk_dir set to NULL")
    chunk_dir <- NULL
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
    
    
    # find the unique combinations of issue_time and horizon at per fold across all the locations
    find_nsamp <- list()
    for(i in names(sigma_kf)){
      find_nsamp[[i]] <- unique(do.call(rbind,unname(lapply(control,function(x){data.table(issue_ind=x$issue_ind[x$kfold==i],horiz_ind=x$horiz_ind[x$kfold==i])}))))
      find_nsamp[[i]] <- find_nsamp[[i]][order(find_nsamp[[i]]$issue_ind, find_nsamp[[i]]$horiz_ind),]
    }
    
    # extract samples calling etr_kf_spatsamp for each fold
    samps <- list() 
    for(i in names(sigma_kf)){
      
      samps[[i]] <- extr_kf_spatsamp(kf_samp_df=find_nsamp[[i]],
                                     uni_kfold = i,
                                     sigmas = sigma_kf,
                                     means = mean_kf,
                                     mvnfast_c = mvnfast_cores,
                                     nsamps = no_samps,
                                     margs = marginals,
                                     chunk = chunk_dir)
    }
    
    
    rm(find_nsamp)
    
    
    
  } else{ if (copulatype=="temporal"){
    
    
    # find number of unique issue_times and lead times per fold across all the locations (use data.table to avoid losing posixct class if present)
    # note that at the moment each location is assumed to have an identical number of lead times - class on cvms!
    find_issue <- list()
    find_horizon <- list()
    for(i in names(sigma_kf)){
      find_issue[[i]] <- unique(do.call(rbind,unname(lapply(control,function(x){data.table(issue_ind=unique(x$issue_ind[x$kfold==i]))}))))
      find_horizon[[i]] <- unique(do.call(rbind,unname(lapply(control,function(x){data.table(horiz_ind=unique(x$horiz_ind[x$kfold==i]))}))))
      find_issue[[i]] <- setorder(find_issue[[i]],issue_ind)
      find_horizon[[i]] <- setorder(find_horizon[[i]],horiz_ind)
    }
    
    # extract samples calling extr_kf_temposamp for each fold
    # output in format kfold$<dt(locs)>
    # note that for chunked calculation samps[[i]] is a dt of location, issuetime, and lead time indexes only
    samps <- list() 
    for(i in names(sigma_kf)){
      
      samps[[i]] <- extr_kf_temposamp(issuetimes = find_issue[[i]],
                                      uni_kfold = i,
                                      leadtimes = find_horizon[[i]],
                                      sigmas = sigma_kf,
                                      means = mean_kf,
                                      mvnfast_c = mvnfast_cores,
                                      nsamps = no_samps,
                                      margs = marginals,
                                      chunk = chunk_dir)
    }
    
    
    rm(find_issue,find_horizon)
    
    
    
  }else{stop("copula type mis-specified")}}
  
  # clear deleted tables
  # not really needed but keeps the memory down before a parallel process
  invisible(gc())
  
  # set control data.tables
  cont_ids <- lapply(control,function(x){data.table(issue_ind=x$issue_ind,horiz_ind=x$horiz_ind,sort_ind=1:length(x$issue_ind))})
  
  print(paste0("Transforming samples into original domain"))
  
  
  if (class(marginals[[1]])[1]%in%c("MultiQR")){
    
    
    method_list <- lapply(control,function(x){x$PIT_method})
    CDFtail_list <- lapply(control,function(x){x$CDFtails})
    
    
    if(!is.null(chunk_dir)){ ### chunked calc - meh
      
      ## get list of locations
      loc_list <- as.list(names(marginals))
      names(loc_list) <- names(marginals)
      
      ## for each kfold inverse pit for each location
      samps_temp <- list()
      for(i in names(sigma_kf)){
        
        print(i)
        
        samps_temp[[i]] <-  mcmapply(function(name_loc, samps_indx,...){
          
          ## get indexes for fold i and location name_loc
          indx <- range(samps_indx[loc_id==name_loc,which=TRUE])
          
          # apply inv pit function loading in the data for fold i and location name_loc
          fastinvpit(dt_samps = fst::read_fst(paste0(chunk_dir,i,"_temp.fst"),from = indx[1], to = indx[2], as.data.table = TRUE),
                     no_samps = no_samps,
                     ...)
          
        },name_loc = loc_list, dt_contr = cont_ids, qrdata = marginals, method = method_list, tails = CDFtail_list,
        MoreArgs = list(samps_indx = samps[[i]],...), mc.cores = mcmapply_cores, SIMPLIFY = FALSE)
        
        ## clear memory before next parallel loop
        invisible(gc())
        
      }
      
      
      for(i in names(sigma_kf)){
        file.remove(paste0(chunk_dir,i,"_temp.fst"))
      }
      
      ## order data.tables properly - from kfold$location to location dts
      samps <- list()
      invisible(gc())
      for(i in names(marginals)){
        
        samps[[i]] <- rbindlist(lapply(samps_temp,function(z){z[[i]][,loc_id:=NULL]}))
        
        ## delete samps temp dt -- all folds location i
        lapply(samps_temp,function(z){delete_memsafe(z[[i]],z[[i]][!is.na(V1),which = TRUE])})
        invisible(gc())
        
      }
      
      rm(samps_temp)
      
      
    } else{ #### no chunked calc
      
      
      samps <- split(rbindlist(samps), by="loc_id",keep.by = FALSE)
      invisible(gc())
      # apply function and transfrom to original domain
      samps <- mcmapply(function(...){fastinvpit(no_samps = no_samps,...)},
                        dt_samps = samps, dt_contr = cont_ids, qrdata = marginals,method = method_list,tails = CDFtail_list,
                        SIMPLIFY = FALSE,mc.cores = mcmapply_cores,MoreArgs = list(...))
      
      
      
    }
    
    
    ## setorder and remove time indx
    lapply(samps,function(x){setorder(x,sort_ind)})
    lapply(samps,function(x){x[,c("sort_ind","issue_ind","horiz_ind"):=NULL]})
    
    
  } else {
    
    ## change to similar to multi QR
    ## add S3 support for PPD?
    samps <- split(rbindlist(samps), by="loc_id",keep.by = FALSE)
    samps <- mapply(merge.data.table,x = cont_ids,y = samps,MoreArgs = list(all.x=T),SIMPLIFY = F)
    rm(cont_ids)
    invisible(gc())
    
    ## preserve order of merged scenario table with input control table
    samps <- lapply(samps,function(x){setorder(x,sort_ind)})
    
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





##### helper functions....

# This function is for extracting temporal/spatio-temporal scenario samples for a single fold
extr_kf_temposamp <- function(issuetimes, uni_kfold, leadtimes, sigmas, means, mvnfast_c, nsamps, margs, chunk){
  
  print(paste0("taking samples for --- ",uni_kfold))
  
  ## take cholesky decomp. so we only have to do it once for all issuetimes in fold
  chol_mat <- chol(sigmas[[uni_kfold]])
  
  arg_list <- list(n=nsamps,sigma=chol_mat,mu=means[[uni_kfold]], isChol = TRUE, ncores = mvnfast_c)
  
  # sample from multivariate gaussian --> list of matrices for each issuetime -->  samples to cols and horizon to rows --> convert to uniform domain --> data.table class --> rbindlist
  kf_samps <- rbindlist(replicate(n=length(issuetimes$issue_ind),expr=data.table(pnorm(t(do.call(eval(parse(text="rmvn")),args=arg_list)))),
                                  simplify = F),
                        idcol = TRUE)
  
  ## add location id column
  kf_samps[,loc_id := rep(names(margs),each = .N/length(margs)),by=.(.id)]
  
  ## rename and redefine .id column to corresponding issuetime
  kf_samps[,.id:=issuetimes$issue_ind[.id]]
  setnames(kf_samps,".id","issue_ind")
  
  ## add leadtime id column 
  kf_samps[,horiz_ind:=leadtimes$horiz_ind,by=.(issue_ind,loc_id)]
  
  # not really needed but good for inspecting the table when coding this :p
  setcolorder(kf_samps,c("issue_ind","loc_id","horiz_ind"))
  
  if(!is.null(chunk)){
    
    ## set order so we can load subsets of bmu_ids via fst
    setorder(kf_samps,loc_id,issue_ind,horiz_ind)
    
    fst::write_fst(kf_samps, path = paste0(chunk,uni_kfold,"_temp.fst"), compress = 100)
    
    kf_samps[,c(paste0("V",1:nsamps)):=NULL]
    invisible(gc())
    
  }
  
  
  return(kf_samps)
  
}



# This function is for extracting spatial scenario samples per fold
extr_kf_spatsamp <- function(kf_samp_df, uni_kfold, sigmas, means, mvnfast_c, nsamps, margs, chunk){
  
  print(paste0("taking samples for --- ",uni_kfold))
  
  ## take cholesky decomp. so we only have to do it once for all issuetimes in fold
  chol_mat <- chol(sigmas[[uni_kfold]])
  
  arg_list <- list(n=nsamps,sigma=chol_mat,mu=means[[uni_kfold]], isChol = TRUE, ncores = mvnfast_c)
  
  # sample from multivariate gaussian --> list of matrices for each row 
  kf_samps <- replicate(n=nrow(kf_samp_df),expr=do.call(eval(parse(text="rmvn")),args=arg_list),simplify = F)
  
  # transform sample rows ---> samples in cols and time_ind/location in rows --> convert to uniform domain
  kf_samps <- lapply(kf_samps,function(x){pnorm(t(x))})
  
  # much more computationally effective to wrap data.table here, because we're taking samples at each row of kf_samp_df
  kf_samps <- data.table(docall("rbind",kf_samps))
  
  kf_samps[,loc_id := rep(names(margs),nrow(kf_samp_df))]
  setkey(kf_samps,loc_id)
  
  # bind with kf_samp_df time indices
  kf_samps <- cbind(rbindlist(replicate(length(margs),kf_samp_df,simplify=FALSE)),kf_samps)
  
  if(!is.null(chunk)){
    
    ## set order so we can load subsets of bmu_ids via fst
    setorder(kf_samps,loc_id,issue_ind,horiz_ind)
    
    fst::write_fst(kf_samps, path = paste0(chunk,uni_kfold,"_temp.fst"), compress = 100)
    
    kf_samps[,c(paste0("V",1:nsamps)):=NULL]
    invisible(gc())
    
  }
  
  
  return(kf_samps)
  
}


### memory safe function for deleting rows due to potential size of datasets
### credit https://stackoverflow.com/questions/10790204/how-to-delete-a-row-by-reference-in-data-table
delete_memsafe <- function(DT, del.idxs) { ## del.idxs here is the rows to remove...
  keep.idxs <- setdiff(DT[, .I], del.idxs);
  cols = names(DT);
  DT.subset <- data.table(DT[[1]][keep.idxs]);
  setnames(DT.subset, cols[1]);
  for (col in cols[2:length(cols)]) {
    DT.subset[, (col) := DT[[col]][keep.idxs]];
    DT[, (col) := NULL];  # delete
  }
  return(DT.subset);
}



# fast inversse pit function
fastinvpit <- function(dt_samps,dt_contr,qrdata, no_samps, ...){
  
  # subsetting by reference -- join on control issue/leadtime --- add sort_ind var from cont_ids
  dt_samps[dt_contr,sort_ind:=sort_ind,on = .(issue_ind,horiz_ind)]
  # delete rows where there is nomatch...by reference someday https://github.com/Rdatatable/data.table/issues/635
  dt_samps <- delete_memsafe(dt_samps, dt_samps[,which(is.na(sort_ind))])
  # clear deleted tables
  invisible(gc())
  ## inverse CDF function for each unique row --> sort by sort_ind
  dt_samps[,paste0("V",1:no_samps):=as.list(contCDF(quantiles = qrdata[sort_ind,],inverse = TRUE,...)(as.numeric(.SD))),
           keyby=.(sort_ind),.SDcols=paste0("V",1:no_samps)]
  
}




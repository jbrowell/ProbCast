#' Generate multivariate forecasts
#'
#' This function produces a list of multivariate scenario forecasts in the marginal domain from the spatial/tempral/spatiotemporal covariance matrices and marginal distributions
#' @param copulatype Either "spatial" of "temporal", note that spatio-temporal can be generated via "temporal" setting
#' @param no_samps Number of scenarios required
#' @details Details go here...
#' @return A list or data frame of multivariate scenario forecasts
#' @export
samps_to_scens <- function(copulatype,no_samps,list_margins,list_sigma,list_mean,control,...){
  
  
  if(length(list_margins)>1){
    if(all(sapply(lapply(list_margins,nrow), identical, lapply(list_margins,nrow)[[1]]))==FALSE){
      stop("margins not of equal size")
    }
  }
  
  if(length(list_margins)>1){
    if(all(sapply(lapply(list_margins,class), identical, lapply(list_margins,class)[[1]]))==FALSE){
      stop("multiple margins must be same class for now...")
    }
  }
  
  
  
  if(copulatype=="spatial"){
    
    samps <- list()
    #### loop over no_of folds
    f_ind <- 1
    for (i in (unique(control$kfold))){
      
      print(i)
      arg_list <- list(n=no_samps,sigma=list_sigma[[f_ind]],mean=list_mean[[f_ind]],...)
      
      #### take no of samples corresponding to no_row in margin
      kf_samps <- replicate(n=sum(control$kfold==i),expr=t(apply(do.call(eval(parse(text="mvtnorm::rmvnorm")),args=arg_list),1,pnorm)))
      
      kf_samps2 <- array(NA,dim = c(dim(kf_samps)[3],dim(kf_samps)[1],length(list_margins)))
      for (j in 1:length(list_margins)){
        
        temp <- kf_samps[,j,]
        temp<- t(temp)
        kf_samps2[,,j] <- temp
        rm(temp)
      }
      
      rm(kf_samps)
      kf_samps <- kf_samps2
      rm(kf_samps2)
      
      
      samps[[i]] <- kf_samps
      rm(kf_samps)
      f_ind <- f_ind+1
    }
  } else{ if (copulatype=="temporal"){
    
    
    samps <- list()
    #### loop over no_of folds
    f_ind <- 1
    for (i in (unique(control$kfold))){
      
      print(i)
      arg_list <- list(n=no_samps,sigma=list_sigma[[f_ind]],mean=list_mean[[f_ind]],...)
      
      #### take no of samples corresponding to no_row in margin
      
      kf_samps <- replicate(n=length(unique(control$issue_ind[control$kfold==i])),
                            expr=t(apply(do.call(eval(parse(text="mvtnorm::rmvnorm")),args=arg_list),1,pnorm)))
      
      if(length(list_margins)>1){
        kf_samps2 <- array(NA,dim = c(length(unique(control$issue_ind[control$kfold==i]))*length(unique(control$horiz_ind)),dim(kf_samps)[1],length(list_margins)))
        kf_ind <- seq(0,ncol(list_sigma[[f_ind]]),length(unique(control$horiz_ind)))
        for (j in 1:length(list_margins)){
          
          temp <- kf_samps[,(kf_ind[j]+1):kf_ind[j+1],]
          dim(temp) <- c(dim(temp)[1],dim(temp)[2]*dim(temp)[3])
          temp<- t(temp)
          kf_samps2[,,j] <- temp
          rm(temp)
        }
        
        rm(kf_samps)
        kf_samps<- kf_samps2
        rm(kf_samps2)
        
      } else{
        
        dim(kf_samps) <- c(dim(kf_samps)[1],dim(kf_samps)[2]*dim(kf_samps)[3])
        kf_samps <- t(kf_samps)
        
      }
      
      
      samps[[i]] <- kf_samps
      rm(kf_samps)
      f_ind <- f_ind+1
    }
    
    
    ###now remove rows which are not present in marginals....
    for (i in (unique(control$kfold))){
      
      mask_hind <- rep(1:(nrow(samps[[i]])/length(unique(control$issue_ind[control$kfold==i]))),length(unique(control$issue_ind[control$kfold==i])))
      mask_iind <- rep(1:length(unique(control$issue_ind[control$kfold==i])),length(unique(control$horiz_ind)))
      mask_iind <- mask_iind[order(mask_iind)]
      
      maskdf <- as.data.frame(cbind(mask_iind,mask_hind))
      
      hind_m <- unique(control$horiz_ind)
      hind_m <- hind_m[order(hind_m)]
      maskdf$mask_hind <- hind_m[maskdf$mask_hind]
      
      iind <- control$issue_ind[control$kfold==i]
      iind_m <- diff(c(0,which(diff(iind)!=0),length(iind)))
      iind_u <- NULL
      for (j in 1:length(unique(control$issue_ind[control$kfold==i]))){
        iind_u <- c(iind_u,rep(j,iind_m[j]))
      }
      
      
      hind <- control$horiz_ind[control$kfold==i]
      maskdf2 <- as.data.frame(cbind(mask_iind=iind_u,mask_hind=hind))
      maskdf2$indy <- 1
      
      
      maskdf <- merge.data.frame(maskdf,maskdf2,all.x = T)
      maskdf <- maskdf[order(maskdf[,1], maskdf[,2]), ]
      
      mask_ind <- which(is.na(maskdf$indy))
      
      if (length(list_margins)>1){
        samps[[i]] <- samps[[i]][-c(mask_ind),,]
      } else{samps[[i]] <- samps[[i]][-c(mask_ind),]}
      rm(mask_ind,mask_hind,mask_iind,maskdf,maskdf2,iind,iind_m,iind_u,hind,hind_m)
      
      
    }
    
    
  }else{stop("copula type mis-specified")}}
  
  
  if (length(dim(samps[[1]]))>2){
    sampstemp <- list()
    for (i in 1:length(list_margins)){
      sampstemp[[i]] <- as.data.frame(matrix(NA,nrow = length(control$kfold),ncol = no_samps))
      for (j in (unique(control$kfold))){
        
        sampstemp[[i]][control$kfold==j,] <- samps[[j]][,,i]
      }
    }
  } else{
    sampstemp <- as.data.frame(matrix(NA,nrow = length(control$kfold),ncol = no_samps))
    for (j in (unique(control$kfold))){
      
      sampstemp[control$kfold==j,] <- samps[[j]]
    }
    
  }
  
  ### transform Unifrom Variable into original domain
  ### add support for PPD
  if (class(list_margins[[1]])[1]%in%c("MultiQR")){
    
    if (length(dim(samps[[1]]))>2){
      sampsfinal <- list()
      for (j in 1:length(list_margins)){
        print(paste0("Transforming samples into original domain --- margin ",j))
        
        sampsfinal[[j]] <- as.data.frame(PIT.MultiQR(qrdata=list_margins[[j]],obs=sampstemp[[j]],inverse=TRUE,method=control$PIT_method,tails=control$CDFtails))
      }
    } else{
      
      print(paste0("Transforming samples into original domain"))
      
      sampsfinal <- as.data.frame(PIT.MultiQR(qrdata=list_margins[[1]],obs=sampstemp,inverse=TRUE,method=control$PIT_method,tails=control$CDFtails))
      
    }
  } else {
    
    
    if (length(dim(samps[[1]]))>2){
      sampsfinal <- list()
      for (j in 1:length(list_margins)){
        print(paste0("Transforming samples into original domain --- margin ",j))
        gooddata <- rowSums(is.na(list_margins[[j]]))==0
        good_params <- list_margins[[j]][gooddata,]
        good_samps <- as.data.frame(lapply(sampstemp[[j]][gooddata,],function(x){do.call(control$q_fun,as.list(cbind(p=x,good_params)))}))
        sampsfinal[[j]] <- as.data.frame(matrix(NA,nrow(list_margins[[j]]),ncol = no_samps))
        sampsfinal[[j]][gooddata,] <- good_samps
      }
    } else{
      
      print(paste0("Transforming samples into original domain"))
      
      gooddata <- rowSums(is.na(list_margins[[1]]))==0
      good_params <- list_margins[[1]][gooddata,]
      good_samps <- as.data.frame(lapply(sampstemp[gooddata,],function(x){do.call(control$q_fun,as.list(cbind(p=x,good_params)))}))
      sampsfinal <- as.data.frame(matrix(NA,nrow(list_margins[[1]]),ncol = no_samps))
      sampsfinal[gooddata,] <- good_samps
    }
    
    
  }
  
  return(sampsfinal)
  
  
}
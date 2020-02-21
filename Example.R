### StrathCast Extended Example
require(devtools)
require(roxygen2)
require(rstudioapi)

PackagePath <- dirname(getActiveDocumentContext()$path)
setwd(PackagePath)

# Update package documentation
document(pkg = ".")
# Install from local repository
install(PackagePath)
# Load Package
require(ProbCast)



setwd(PackagePath)
### Testing Functionalisty of StrathCast #####

## Add some features first...

Wind$WS100 <- sqrt(Wind$U100^2+Wind$V100^2)
Wind$Power <- pmin(Wind$WS100,11)^3

## Set-up simple kfold CV
Wind$kfold <- c(rep(c("fold1","fold2","fold3"),each=4000),rep("Test",nrow(Wind)-12000))
# Wind$kfold <- c(rep(c("fold1","fold2","fold3","fold4"),each=nrow(Wind)/4))

### Multiple Quantile Regression using GBM ####
test1<-list(data=Wind)

# test1$gbm_mqr <- MQR_gbm(data = test1$data,
#                          formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
#                          gbm_params = list(interaction.depth = 3,
#                                            n.trees = 1000,
#                                            shrinkage = 0.05,
#                                            n.minobsinnode = 20,
#                                            bag.fraction = 0.5,
#                                            keep.data = F),
#                          quantiles = seq(0.1,0.9,by=0.1),
#                          Sort = T,
#                          SortLimits = list(U=0.999,L=0.001),
#                          pred_ntree = 1000)
# 
# test1$gbm_mqr <- MQR_gbm(data = test1$data,
#                          formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
#                          gbm_params = list(interaction.depth = 3,
#                                            n.trees = 1000,
#                                            shrinkage = 0.05,
#                                            n.minobsinnode = 20,
#                                            bag.fraction = 0.5,
#                                            keep.data = F),
#                          parallel = T,
#                          cores = 3,
#                          quantiles = seq(0.1,0.9,by=0.1),
#                          Sort = T,
#                          SortLimits = list(U=0.999,L=0.001),
#                          pred_ntree = 1000)


test1$gbm_mqr <- MQR_gbm(data = test1$data,
                         formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                         gbm_params = list(interaction.depth = 3,
                                           n.trees = 1000,
                                           shrinkage = 0.05,
                                           n.minobsinnode = 20,
                                           bag.fraction = 0.5,
                                           keep.data = F),
                         parallel = T,
                         cores = detectCores(),
                         quantiles = seq(0.05,0.95,by=0.05),
                         Sort = T,
                         SortLimits = list(U=1,L=0),
                         pred_ntree = 1000,
                         para_over_q = T)


# test1$gbm_mqr <- test1$gbm_mqr3

par(mar=c(3,3,0.5,1))  # Trim margin around plot [b,l,t,r]
par(tcl=0.35)  # Switch tick marks to insides of axes
par(mgp=c(1.5,0.2,0))  # Set margin lines; default c(3,1,0) [title,labels,line]
par(xaxs="r",yaxs="r")  # Extend axis limits by 4% ("i" does no extension)

# plot(test1$gbm_mqr[1:72,],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1); axis(1,0:12*6,pos=-0.07); axis(2,las=1)
# lines(test1$data$TARGETVAR[1:250],lwd=3)

## thought it might be good to show 1 issueTime? obv a good one :p
# set.seed(222)
# i_ts <- sample(unique(test1$data$ISSUEdtm),1)

i_ts <- unique(test1$data$ISSUEdtm)[3]

plot(test1$gbm_mqr[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
# lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

reliability(qrdata = test1$gbm_mqr,
            realisations = test1$data$TARGETVAR,
            kfolds = test1$data$kfold)

pinball(qrdata = test1$gbm_mqr,
         realisations = test1$data$TARGETVAR,
         kfolds = test1$data$kfold)


reliability(qrdata = test1$gbm_mqr,
            realisations = test1$data$TARGETVAR)

reliability(qrdata = test1$gbm_mqr[test1$data$kfold=="Test",],
            realisations = test1$data$TARGETVAR[test1$data$kfold=="Test"],
            subsets = test1$data$WS100[test1$data$kfold=="Test"],
            breaks = 4,
            bootstrap = 100)

pinball(qrdata = test1$gbm_mqr,
        realisations = test1$data$TARGETVAR)

index <- 54
x <- seq(0,1,by=0.001)
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],method = "spline")
plot(x,cdf(x),type="l",xlab="Target Variable",ylab="CDF",axes=F); axis(1); axis(2,las=2); #grid()
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],method = "linear")
lines(x,cdf(x),lty=2,col=2)
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],kfold = NA,method = "spline", tails=list(method="exponential",L=0,U=1,nBins=5,preds=test1$gbm_mqr,targetvar=test1$data$TARGETVAR,ntailpoints=25))
lines(x,cdf(seq(0,1,by=0.001)),lty=3,col=4)
points(test1$gbm_mqr[index,],as.numeric(gsub("q","",colnames(test1$gbm_mqr[index,])))/100)
legend(0.01,1,c("Predicted Quantiles","Linear","Spline","Spline with Exponential Tails"),
       pch=c(1,NA,NA,NA),lty=c(NA,2,1,3),col=c(1,2,1,4),bty="n")

# test1$X_gbm <- PIT(test1$gbm_mqr,test1$data$TARGETVAR,method = "spline",tails=list(method="exponential",L=0,U=1,nBins=5,preds=test1$gbm_mqr,targetvar=test1$data$TARGETVAR,ntailpoints=25))
test1$X_gbm <- PIT(test1$gbm_mqr,test1$data$TARGETVAR,method = "spline",tails=list(method="interpolate",L=0,U=1))
hist(test1$X_gbm,breaks = 50,freq=F,ylim = c(0,3)); lines(c(0,1),c(1,1),lty=2)

### Parametric PredDist Using GAMLSS ####


test1$ppd <- Para_gamlss(data = test1$data,
                         formula = TARGETVAR~bs(WS100,df=3),
                         sigma.formula = ~WS100,
                         sigma.start = 0.05,
                         nu.formula = ~WS100,
                         tau.formula = ~WS100,
                         family = BEINF, # NO
                         method=mixed(20,10))


summary(test1$ppd$fold1)
plot(test1$ppd$fold1)

test1$gamlssParams <- PPD_2_MultiQR(data=test1$data,
                                   models = test1$ppd,
                                   params = T)

# some issue with the gamlss predictions here, needs futher digging...
test1$gamlssParams[which(test1$gamlssParams[,1]>=1),1] <- 0.99999
test1$gamlssParams[which(test1$gamlssParams[,2]>=1),2] <- 0.99999

test1$gamlss_mqr <- PPD_2_MultiQR(data=test1$data,
                                  models = test1$ppd,
                                  params = F)

plot(test1$gamlss_mqr[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

reliability(qrdata = test1$gamlss_mqr,
            realisations = test1$data$TARGETVAR,
            kfolds = test1$data$kfold)

pinball(qrdata = test1$gamlss_mqr,
         realisations = test1$data$TARGETVAR,
         kfolds = test1$data$kfold)


# test1$data[test1$data$TARGETVAR<0 | test1$data$TARGETVAR>1,]
test1$X_gamlss <- PIT(test1$ppd,data = test1$data)
hist(test1$X_gamlss,breaks = 50,freq = F,ylim = c(0,3)); lines(c(0,1),c(1,1),lty=2)

#######################
#### generate temporal scenarios using the gaussion copula and gbm_MQR marginals
#######################


# define temporal covariance matrix
u_obsind <- data.frame(kfold=test1$data$kfold,lead_time=as.numeric(test1$data$TARGETdtm-test1$data$ISSUEdtm),i_time=test1$data$ISSUEdtm,u_obs = test1$X_gbm)
u_obswide <- reshape(u_obsind,idvar = "i_time",direction = "wide",v.names = "u_obs",timevar = "lead_time",sep = "_")
u_obswide <- u_obswide[order(u_obswide$i_time),]

# function doesn't use "Test" data when defining any of the matrices if kfold is specified
cvm_gbm <- covcor_matrix(u_data = u_obswide[,-c(1:2)],cov_cor = "covariance",kfold = u_obswide$kfold, scale = T, method = "pearson")
#### requires lattice --- hour 24 looks weird --- are you sure the target_time goes from 1-24 and not 0-23 for each issue time?
col6 <- colorRampPalette(c("blue","cyan","yellow","red"))
lattice::levelplot(cvm_gbm[["Test"]], xlab="lead time [hours]", ylab="lead time [hours]",col.regions=col6(600), cuts=100, at=seq(-0.3,1,0.01),
                   scales=list(x=list(at=seq(0,24,3),rot=90),y=list(at=seq(0,24,3)),tck=0.3,cex=1.1),
                   main="Test --- Covariance")

# sample cvm and convert to power domain
f_nsamp <- 200
mean_list <- list()
for (i in levels(unique(u_obsind$kfold))){
  mean_list[[i]] <- rep(0, 24)
}
## method for parametric pred dist.
# scen_gbm <- samps_to_scens(copulatype = "temporal",no_samps = f_nsamp,marginals = list(loc_1 = test1$gbm_mqr),sigma_kf = cvm_gbm,mean_kf = mean_list,
#                           control=list(loc_1 = list(kfold = u_obsind$kfold,issue_ind=u_obsind$i_time,horiz_ind=u_obsind$lead_time,
#                                                     PIT_method="spline",
#                                                     CDFtails = list(method="interpolate",L=0,U=1,ntailpoints=100))))
scen_gbm <- samps_to_scens(copulatype = "temporal",no_samps = f_nsamp,marginals = list(loc_1 = test1$gbm_mqr),sigma_kf = cvm_gbm,mean_kf = mean_list,
                           control=list(loc_1 = list(kfold = u_obsind$kfold,issue_ind=u_obsind$i_time,horiz_ind=u_obsind$lead_time,
                                                     PIT_method="linear",
                                                     CDFtails = list(method="interpolate",L=0,U=1,ntailpoints=100))))



matplot(scen_gbm$loc_1[which(test1$data$ISSUEdtm==i_ts),],type="l",ylim=c(0,1),lty=1,
        xlab="Lead Time [Hours]",ylab="Power [Capacity Factor]",
        col=gray(0.1,alpha = 0.1),axes = F); axis(1,1:24,pos=-0.07); axis(2,las=1)
# lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=2)
# legend("bottomleft",c("scenarios","measured"),col = c("grey75","black"),pch=c(NA,NA,NA),bty="n",lty=1)



#######################
#### generate temporal scenarios using the gaussion copula and PPD marginals
#######################

# define temporal covariance matrix
u_obsind <- data.frame(kfold=test1$data$kfold,lead_time=as.numeric(test1$data$TARGETdtm-test1$data$ISSUEdtm),i_time=test1$data$ISSUEdtm,u_obs = test1$X_gamlss)
u_obswide <- reshape(u_obsind,idvar = "i_time",direction = "wide",v.names = "u_obs",timevar = "lead_time",sep = "_")
u_obswide <- u_obswide[order(u_obswide$i_time),]

# function doesn't use "Test" data when defining any of the matrices if kfold is specified
cvm_gamlss <- covcor_matrix(u_data = u_obswide[,-c(1:2)],cov_cor = "covariance",kfold = u_obswide$kfold, scale = T, method = "pearson")
#### requires lattice --- hour 24 looks weird --- are you sure the target_time goes from 1-24 and not 0-23 for each issue time?
col6 <- colorRampPalette(c("blue","cyan","yellow","red"))
lattice::levelplot(cvm_gamlss[["Test"]], xlab="lead time [hours]", ylab="lead time [hours]",col.regions=col6(600), cuts=100, at=seq(-0.3,1,0.01),
                   scales=list(x=list(at=seq(0,24,3),rot=90),y=list(at=seq(0,24,3)),tck=0.3,cex=1.1),
                   main="Test --- Covariance")

# sample cvm and convert to power domain
# method for parametric pred dist.
scen_gamlss <- samps_to_scens(copulatype = "temporal",no_samps = f_nsamp,marginals = list(loc_1 = test1$gamlssParams),sigma_kf = cvm_gamlss,mean_kf = mean_list,
                           control=list(loc_1 = list(kfold = u_obsind$kfold,issue_ind=u_obsind$i_time,horiz_ind=u_obsind$lead_time,
                                                     q_fun = gamlss.dist::qBEINF)))

matplot(scen_gamlss$loc_1[which(test1$data$ISSUEdtm==i_ts),],type="l",ylim=c(0,1),lty=1,
        xlab="Lead Time [Hours]",ylab="Power [Capacity Factor]",
        col=gray(0.1,alpha = 0.1),axes = F); axis(1,1:24,pos=-0.07); axis(2,las=1)
# lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=2)
# legend("topleft",c("scenarios","measured"),col = c("grey75","black"),pch=c(NA,NA,NA),bty="n",lty=1)





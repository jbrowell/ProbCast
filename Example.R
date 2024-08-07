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



### Testing Functionality of ProbCast #####

## Add some features first...

Wind$WS100 <- sqrt(Wind$U100^2+Wind$V100^2)
Wind$WD100 <- atan2(Wind$V100,Wind$U100)
Wind$Power <- pmin(Wind$WS100,9.5)^3 / 9.5^3

## Set-up simple kfold CV. NB --- For scenario forecasting make sure the CV folds don't cross issue times
Wind$kfold <- "Fold 1"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2012-06-30",tz="UTC")] <- "Fold 2"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2012-12-31",tz="UTC")] <- "Fold 3"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2013-06-30",tz="UTC")] <- "Test"


### Multiple linear quantile regression with MQR_rq
# require(splines2)
model_rq = qreg_mrq(data=Wind,
                    formula = TARGETVAR~1+Power+WS100,
                    # formula = TARGETVAR~1+bSpline(WS100,df=10) +
                    #   bSpline(WS100,df=10) + 
                    #   bSpline(WD100,df = 6,periodic = T,Boundary.knots = c(-pi,pi)),
                    quantiles = 1:19/20,
                    offset = "Power",
                    cv_folds = "kfold",
                    sort_limits = list(U=1,L=0))

plot(model_rq$mqr_pred[1:100+sample(1:16000,1),])

reliability(qrdata = model_rq$mqr_pred,
            realisations = Wind$TARGETVAR,
            kfolds = Wind$kfold)

pinball(qrdata = model_rq$mqr_pred,
            realisations = Wind$TARGETVAR,
            kfolds = Wind$kfold)


plot(Width~Interval,sharpness(qrdata = model_rq$mqr_pred,
                              realisations = Wind$TARGETVAR),
     type="b")



### Multiple Quantile Regression using lightGBM ####
lgbm_model <- qreg_lightgbm(data=Wind,
                            formula=TARGETVAR~U100+V100+U10+V10+WS100,
                            quantiles = seq(0.05,0.95,by=0.05),
                            cv_folds = "kfold",
                            sort = TRUE,
                            sort_limits = list(L=0, U=1),
                            cores=detectCores() - 1,
                            lightgbm_params = list(n_estimators=1000,
                                                   learning_rate=0.01,
                                                   min_data_in_leaf=30,
                                                   bagging_freq=1,
                                                   bagging_fraction=0.9,
                                                   min_gain_to_split=0.25,
                                                   lambda_l1=10,
                                                   lambda_l2=10,
                                                   force_col_wise=TRUE))


### re-train models over the Test set, every 28 days
updated_lgbm_model <- retrain_all(lgbm_model,
                                  data=Wind,
                                  retrain_daily_frequency=28,
                                  issue_datetime_column='ISSUEdtm',
                                  cv_folds='kfold',
                                  cores=detectCores() - 1)

plot_idx = 1:100+sample(1:16000,1)
plot(updated_lgbm_model$mqr_pred[plot_idx,])
lines(1:100, Wind[plot_idx, 'TARGETVAR'])

reliability(qrdata = updated_lgbm_model$mqr_pred,
            realisations = Wind$TARGETVAR,
            kfolds = Wind$kfold)

pinball(qrdata = updated_lgbm_model$mqr_pred,
        realisations = Wind$TARGETVAR,
        kfolds = Wind$kfold)


plot(Width~Interval,sharpness(qrdata = updated_lgbm_model$mqr_pred,
                              realisations = Wind$TARGETVAR),
     type="b")


### Multiple Quantile Regression using GBM ####
test1<-list(data=Wind)
test1$gbm_mqr <- qreg_gbm(data = test1$data,
                          formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                          cv_folds = "kfold",
                          interaction.depth = 3,
                          n.trees = 1000,
                          shrinkage = 0.05,
                          n.minobsinnode = 20,
                          bag.fraction = .5,
                          keep.data = F,
                          quantiles = seq(0.05,0.95,by=0.05),
                          sort = T,
                          sort_limits = list(U=1,L=0),
                          pred_ntree = 1000,
                          cores=detectCores(),
                          only_mqr = TRUE)



par(mar=c(3,3,0.5,1))  # Trim margin around plot [b,l,t,r]
par(tcl=0.35)  # Switch tick marks to insides of axes
par(mgp=c(1.5,0.2,0))  # Set margin lines; default c(3,1,0) [title,labels,line]
par(xaxs="r",yaxs="r")  # Extend axis limits by 4% ("i" does no extension)

i_ts <- unique(test1$data$ISSUEdtm)[3]

plot(test1$gbm_mqr[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

reliability(qrdata = test1$gbm_mqr,
            realisations = test1$data$TARGETVAR,
            subsets = test1$data$WS100,
            breaks = c(4,10),bootstrap = 100)

pinball(qrdata = test1$gbm_mqr,
        realisations = test1$data$TARGETVAR,
        kfolds = test1$data$kfold)


reliability(qrdata = test1$gbm_mqr,
            realisations = test1$data$TARGETVAR,bootstrap = 100)

reliability(qrdata = test1$gbm_mqr[test1$data$kfold=="Test",],
            realisations = test1$data$TARGETVAR[test1$data$kfold=="Test"],
            subsets = test1$data$WS100[test1$data$kfold=="Test"],
            breaks = 4,
            bootstrap = 100)

pinball(qrdata = test1$gbm_mqr,
        realisations = test1$data$TARGETVAR,
        bootstrap = 100)

pinball(qrdata = test1$gbm_mqr,
        realisations = test1$data$TARGETVAR,
        kfolds = test1$data$kfold,
        bootstrap = 100,ylim=c(0,.08))

pinball(qrdata = test1$gbm_mqr[test1$data$kfold=="Test",],
        realisations = test1$data$TARGETVAR[test1$data$kfold=="Test"],
        subsets = test1$data$WS100[test1$data$kfold=="Test"],
        breaks = 4,
        bootstrap = 100,
        ylim=c(0,.1))

pinball(qrdata = test1$gbm_mqr[test1$data$kfold=="Test",],
        realisations = test1$data$TARGETVAR[test1$data$kfold=="Test"],
        subsets = as.factor((test1$data$TARGETdtm-test1$data$ISSUEdtm)[test1$data$kfold=="Test"]),
        ylim=c(0,.1))





index <- 54
x <- seq(0,1,by=0.001)
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],method = "spline")
plot(x,cdf(x),type="l",xlab="Target Variable",ylab="CDF",axes=F); axis(1); axis(2,las=2); #grid()
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],method = "linear")
lines(x,cdf(x),lty=2,col=2)
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],method = "spline", tails=list(method="dyn_exponential",ntailpoints=25))
lines(x,cdf(x),lty=4,col=5)
points(test1$gbm_mqr[index,],as.numeric(gsub("q","",colnames(test1$gbm_mqr[index,])))/100)
legend(0.01,1,c("Predicted Quantiles","Linear","Spline","Spline with Exponential Tails"),
       pch=c(1,NA,NA,NA),lty=c(NA,2,1,3),col=c(1,2,1,4),bty="n")

## Check inverse matches:
x <- seq(0,1,by=0.001)
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],method = "spline")
plot(x,cdf(x),type="l",xlab="Target Variable",ylab="CDF",axes=F); axis(1); axis(2,las=2); #grid()
inv_cdf <- contCDF(quantiles = test1$gbm_mqr[index,],method = "spline",inverse = T)
lines(inv_cdf(x),x,lty=2,col=2)



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
                         family =  BEINF, #NO,  #
                         method=mixed(20,10))


summary(test1$ppd$`Fold 1`)
plot(test1$ppd$`Fold 1`)

test1$gamlssParams <- PPD_2_MultiQR(data=test1$data,
                                    models = test1$ppd,
                                    params = T)




# some issue with the gamlss predictions here, needs futher digging...
test1$gamlssParams[which(test1$gamlssParams[,1]>=1),1] <- 0.99999
test1$gamlssParams[which(test1$gamlssParams[,2]>=1),2] <- 0.99999

test1$gamlss_mqr <- PPD_2_MultiQR(data=test1$data,
                                  models = test1$ppd,
                                  params = F)

plot(test1$gamlss_mqr[which(test1$data$ISSUEdtm==i_ts),],xlab="Lead time [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
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
for (i in unique(u_obsind$kfold)){
  mean_list[[i]] <- rep(0, 24)
}


## method for gbm pred dist.
scen_gbm <- samps_to_scens(copulatype = "temporal",no_samps = f_nsamp,marginals = list(loc_1 = test1$gbm_mqr),sigma_kf = cvm_gbm,mean_kf = mean_list,
                           control=list(loc_1 = list(kfold = u_obsind$kfold,issue_ind=u_obsind$i_time,horiz_ind=u_obsind$lead_time,
                                                     PIT_method="spline",
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



#######################
#### Evaluate scenarios forecasts using scoringRules
#######################

library(scoringRules)
library(data.table)

# weight matrix function for variogram score
mat <- function(d,horizon){w_vs <- matrix(NA, nrow = d, ncol = d)
for(d1 in 1:d){for(d2 in 1:d){w_vs[d1,d2] <- 0.5^abs(horizon[d1]-horizon[d2])}}
return(w_vs)}



### gbm
FCs <- data.table(cbind(test1$data,scen_gbm$loc_1))
FCs[,horiz:=as.numeric(TARGETdtm - ISSUEdtm)]
test1$mvscore_gbm <- FCs[,list(ES=es_sample(y=TARGETVAR,dat=(as.matrix(.SD))),
                               wVS1=vs_sample(y=TARGETVAR,dat=(as.matrix(.SD)),w=mat(d = .N,horizon = horiz),p=1),
                               wVS.5=vs_sample(y=TARGETVAR,dat=(as.matrix(.SD)),w=mat(d = .N,horizon = horiz),p=.5))
                         ,.SDcols=paste0("scen_",1:f_nsamp),by=c("kfold","ISSUEdtm")]


### gamlss
FCs <- data.table(cbind(test1$data,scen_gamlss$loc_1))
FCs[,horiz:=as.numeric(TARGETdtm - ISSUEdtm)]
test1$mvscore_gamlss <- FCs[,list(ES=es_sample(y=TARGETVAR,dat=(as.matrix(.SD))),
                                  wVS1=vs_sample(y=TARGETVAR,dat=(as.matrix(.SD)),w=mat(d = .N,horizon = horiz),p=1),
                                  wVS.5=vs_sample(y=TARGETVAR,dat=(as.matrix(.SD)),w=mat(d = .N,horizon = horiz),p=.5))
                            ,.SDcols=paste0("scen_",1:f_nsamp),by=c("kfold","ISSUEdtm")]



test1$mvscore_gbm[,lapply(.SD,function(x){mean(x,na.rm = T)}),.SDcols=c("ES","wVS1","wVS.5"),by=.(kfold)]
test1$mvscore_gamlss[,lapply(.SD,function(x){mean(x,na.rm = T)}),.SDcols=c("ES","wVS1","wVS.5"),by=.(kfold)]



# Block bootstrap sampling --- accounting for the temporal correlation of weather patters. Blocks of 7 days...

test1$mvscore_gbm[,block:=as.numeric(floor((ISSUEdtm-as.POSIXct("2012-01-01 00:00:00",tz="UTC"))/(60*60*24*7)))]
test1$mvscore_gamlss[,block:=as.numeric(floor((ISSUEdtm-as.POSIXct("2012-01-01 00:00:00",tz="UTC"))/(60*60*24*7)))]

mv_dt <- rbindlist(list(gbm = test1$mvscore_gbm,gamlss = test1$mvscore_gamlss),idcol = "marginal")
mv_dt[,marginal:=factor(marginal,levels = c("gamlss","gbm"))]
setorder(mv_dt,ISSUEdtm)


evalplot_block <- function(data_table, block,nboot = 100, na.rm = TRUE,score = "ES",...) {
  
  boot <- NULL
  for(i in 1:nboot) {
    bootind <- sample(unique(block), replace = TRUE)
    data <- rbindlist(lapply(bootind,function(x){data_table[block==x]}))
    
    boot <- rbind(boot, data[,as.list(colMeans(.SD,na.rm = na.rm)),.SDcols = score,by=.(marginal)])
    rm(data)
  }
  
  boxplot(data = boot, as.formula(paste0(score,"~ marginal")),ylab = score,xlab = "", ...)
}


### ES CV - gbm

par(mfrow = c(1,2), mar = c(1.5,3,0.5,0),tcl=0.35, mgp=c(1.5,0.2,0), xaxs="r",yaxs="r")
evalplot_block(mv_dt[kfold=="Test"],block = mv_dt[kfold=="Test",block],axes=F,ylim=c(0.4,0.75))
axis(2,seq(0.4,.75,0.05),lwd=2, cex=1.2);axis(1, at=1:2,labels = c("gamlss","gbm"),lwd=2, cex=1.2)
evalplot_block(mv_dt[kfold=="Test"],block = mv_dt[kfold=="Test",block],score="wVS1",axes=F,ylim=c(0.4,0.75))
axis(2,seq(0.4,.75,0.05),lwd=2, cex=1.2);axis(1, at=1:2,labels = c("gamlss","gbm"),lwd=2, cex=1.2)




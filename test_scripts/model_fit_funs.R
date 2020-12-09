########################### temporary example file to test functions and ensure updates are ok


  
rm(list=ls())
.rs.restartR()

# library(foreach)
library(devtools)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

load_all("../", TRUE)


Wind <- ProbCast::Wind

Wind$WS100 <- sqrt(Wind$U100^2+Wind$V100^2)
Wind$Power <- pmin(Wind$WS100,11)^3

## Set-up simple kfold CV. NB --- For scenario forecasting make sure the CV folds don't cross issue times
Wind$kfold <- "Fold 1"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2012-06-30",tz="UTC")] <- "Fold 2"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2012-12-31",tz="UTC")] <- "Fold 3"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2013-06-30",tz="UTC")] <- "Test"

#################################################################################
### Multiple Quantile Regression using GBM ####
#################################################################################

test1<-list(data=Wind)

## test in series
test1$gbm_mqr <- qreg_gbm(data = test1$data,
                          formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                          cv_folds = "kfold",
                          interaction.depth = 3,
                          n.trees = 100,
                          shrinkage = 0.05,
                          n.minobsinnode = 20,
                          bag.fraction = 1,
                          keep.data = F,
                          quantiles = seq(0.1,0.9,by=0.1),
                          sort = T,
                          sort_limits = list(U=0.999,L=0.001),
                          pred_ntree = 100)

# test with modified hyperparameters
test1$gbm_mqr2 <- qreg_gbm(data = test1$data,
                           formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                           cv_folds = "kfold",
                           interaction.depth = 5,
                           n.trees = 100,
                           shrinkage = 0.5,
                           n.minobsinnode = 50,
                           bag.fraction = 1,
                           keep.data = F,
                           quantiles = seq(0.1,0.9,by=0.1),
                           sort = T,
                           sort_limits = list(U=0.999,L=0.001),
                           pred_ntree = 100)

# test in parallel ---> should equal model 1
test1$gbm_mqr3 <- qreg_gbm(data = test1$data,
                           formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                           cv_folds = "kfold",
                           interaction.depth = 3,
                           n.trees = 100,
                           shrinkage = 0.05,
                           n.minobsinnode = 20,
                           bag.fraction = 1,
                           keep.data = F,
                           quantiles = seq(0.1,0.9,by=0.1),
                           sort = T,
                           sort_limits = list(U=0.999,L=0.001),
                           pred_ntree = 100,
                           cores=8L)

# test with save model function 1 worker
tmp <- qreg_gbm(data = test1$data,
                formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                cv_folds = "kfold",
                interaction.depth = 3,
                n.trees = 100,
                shrinkage = 0.05,
                n.minobsinnode = 20,
                bag.fraction = 1,
                keep.data = F,
                quantiles = seq(0.1,0.9,by=0.1),
                sort = T,
                sort_limits = list(U=0.999,L=0.001),
                pred_ntree = 100,
                save_models_path = "./tmp")
rm(tmp)

file.remove(list.files(pattern = "*.rda"))

# test with save model function 8 workers
tmp <- qreg_gbm(data = test1$data,
                formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                cv_folds = "kfold",
                interaction.depth = 3,
                n.trees = 100,
                shrinkage = 0.05,
                n.minobsinnode = 20,
                bag.fraction = 1,
                keep.data = F,
                quantiles = seq(0.1,0.9,by=0.1),
                sort = T,
                sort_limits = list(U=0.999,L=0.001),
                pred_ntree = 100,
                save_models_path = "./tmp",
                cores=8L,
                only_mqr = TRUE)

file.remove(list.files(pattern = "*.rda"))
rm(tmp)





#####################

## back compatibility test
test1$gbm_mqr4 <- MQR_gbm(data = test1$data,
                          formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                          gbm_params = list(interaction.depth = 3,
                                            n.trees = 100,
                                            shrinkage = 0.05,
                                            n.minobsinnode = 20,
                                            bag.fraction = 1,
                                            keep.data = F),
                          quantiles = seq(0.1,0.9,by=0.1),
                          Sort = T,
                          SortLimits = list(U=0.999,L=0.001),
                          pred_ntree = 100)



test1$gbm_mqr5 <- MQR_gbm(data = test1$data,
                          formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                          gbm_params = list(interaction.depth = 5,
                                            n.trees = 100,
                                            shrinkage = 0.5,
                                            n.minobsinnode = 50,
                                            bag.fraction = 1,
                                            keep.data = F),
                          quantiles = seq(0.1,0.9,by=0.1),
                          Sort = T,
                          SortLimits = list(U=0.999,L=0.001),
                          pred_ntree = 100)



test1$gbm_mqr6 <- MQR_gbm(data = test1$data,
                          formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                          gbm_params = list(interaction.depth = 3,
                                            n.trees = 100,
                                            shrinkage = 0.05,
                                            n.minobsinnode = 20,
                                            bag.fraction = 1,
                                            keep.data = F),
                          quantiles = seq(0.1,0.9,by=0.1),
                          Sort = T,
                          SortLimits = list(U=0.999,L=0.001),
                          pred_ntree = 100,
                          cores = 8,
                          parallel = TRUE,
                          para_over_q = FALSE)


### default - autocv
test1$gbm_mqr7 <- qreg_gbm(data = test1$data,
                           formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                           quantiles = seq(0.1,0.9,by=0.1),
                           sort = T,
                           sort_limits = list(U=0.999,L=0.001),
                           pred_ntree = 100)

class(test1$gbm_mqr7$qreg_list)
names(test1$gbm_mqr7$qreg_list)


tmp <- predict(test1$gbm_mqr7$qreg_list,
               newdata = test1$data,
               pred_ntree = 100,
               sort = T,
               sort_limits = list(U=0.999,L=0.001))

data.table(tmp)


tmp <- predict(test1$gbm_mqr7$qreg_list,
               newdata = test1$data,
               which_model = "fold1",
               pred_ntree = 100,
               sort = T,
               sort_limits = list(U=0.999,L=0.001))

data.table(tmp)


tmp <- predict(test1$gbm_mqr7$qreg_list,
               newdata = test1$data,
               which_model = "Fold X",
               pred_ntree = 100,
               sort = T,
               sort_limits = list(U=0.999,L=0.001))




##### OPERATIONAL MODEL
test1$gbm_mqr8 <- qreg_gbm(data = test1$data,
                           formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                           cv_folds = NULL,
                           interaction.depth = 3,
                           n.trees = 100,
                           shrinkage = 0.05,
                           n.minobsinnode = 20,
                           bag.fraction = 1,
                           keep.data = F,
                           quantiles = seq(0.1,0.9,by=0.1))


class(test1$gbm_mqr8$qreg_list)


tmp <- predict(test1$gbm_mqr8$qreg_list,
               newdata = test1$data,
               pred_ntree = 100,
               sort = T,
               sort_limits = list(U=0.999,L=0.001))

data.table(tmp)



tmp <- predict(test1$gbm_mqr8$qreg_list,
               newdata = test1$data,
               which_model = "all_data",
               pred_ntree = 100,
               sort = T,
               sort_limits = list(U=0.999,L=0.001))

data.table(tmp)

tmp <- predict(test1$gbm_mqr8$qreg_list,
               newdata = test1$data,
               pred_ntree = 1,
               quantiles = c(.1,.55,.9))


tmp <- predict(test1$gbm_mqr8$qreg_list,
               newdata = test1$data,
               pred_ntree = 1,
               quantiles = c(.1,.5,.9))


data.table(tmp)

length(test1$gbm_mqr8$qreg_list$all_data$q10$fit)==length(na.omit(test1$data$TARGETVAR))



##### auto cv return only mqr
test1$gbm_mqr9 <- qreg_gbm(data = test1$data,
                           formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                           cv_folds = 5, 
                           interaction.depth = 3,
                           n.trees = 100,
                           shrinkage = 0.05,
                           n.minobsinnode = 20,
                           bag.fraction = 1,
                           keep.data = F,
                           quantiles = seq(0.1,0.9,by=0.1),
                           sort = T,
                           sort_limits = list(U=0.999,L=0.001),
                           pred_ntree = 100,
                           cores=8L,
                           only_mqr = TRUE)


plot(test1$gbm_mqr9[1:100,])






par(mar=c(3,3,0.5,1))  # Trim margin around plot [b,l,t,r]
par(tcl=0.35)  # Switch tick marks to insides of axes
par(mgp=c(1.5,0.2,0))  # Set margin lines; default c(3,1,0) [title,labels,line]
par(xaxs="r",yaxs="r")  # Extend axis limits by 4% ("i" does no extension)

i_ts <- unique(test1$data$ISSUEdtm)[3]

plot(test1$gbm_mqr$mqr_kfold[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

plot(test1$gbm_mqr3$mqr_kfold[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

plot(test1$gbm_mqr4[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

plot(test1$gbm_mqr6[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

plot(test1$gbm_mqr7$mqr_kfold[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)



plot(test1$gbm_mqr2$mqr_kfold[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)


plot(test1$gbm_mqr5[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)






rel <- rbindlist(lapply(names(test1)[-1][c(1,3,4,6,7,2,5)],function(x){
  
  if(class(test1[[x]])[1]!="MultiQR"){
    
    reliability(qrdata = test1[[x]]$mqr_kfold,
                realisations = test1$data$TARGETVAR,
                kfolds = test1$data$kfold,main=x)
    
    
  } else{
    
    
    reliability(qrdata = test1[[x]],
                realisations = test1$data$TARGETVAR,
                kfolds = test1$data$kfold,main=x)
    
    
    
  }
  
  
  
}),idcol="model")



pball <- rbindlist(lapply(names(test1)[-1][c(1,3,4,6,7,2,5)],function(x){
  
  if(class(test1[[x]])[1]!="MultiQR"){
  
  pinball(qrdata = test1[[x]]$mqr_kfold,
          realisations = test1$data$TARGETVAR,
          kfolds = test1$data$kfold,main=x)
    
  } else{
    
    pinball(qrdata = test1[[x]],
            realisations = test1$data$TARGETVAR,
            kfolds = test1$data$kfold,main=x)
    
  }
  
  
}),idcol="model")


fwrite(rel,"reltest.csv")

fwrite(pball,"pballtest.csv")

file.remove(list.files(pattern = "*.csv"))

#################################################################################
### Multiple Quantile Regression using QGAM ####
#################################################################################

load_all("./", TRUE)

test1<-list(data=Wind)

## test in series
test1$gam_mqr <- mqr_qreg_mboost(data = test1$data,
                                formula = TARGETVAR~bbs(WS100,knots=8),
                                quantiles = seq(0.1,0.9,by=0.1),
                                sort = T,
                                sort_limits = list(U=0.999,L=0.001))

# test with modified hyperparameters
test1$gam_mqr2 <- mqr_qreg_mboost(data = test1$data,
                                 formula = TARGETVAR~bbs(WS100,knots=8),
                                 control = mboost::boost_control(mstop=100,nu=0.5),
                                 quantiles = seq(0.1,0.9,by=0.1),
                                 sort = T,
                                 sort_limits = list(U=0.999,L=0.001))


# test in parallel ---> should equal model 1
test1$gam_mqr3 <- mqr_qreg_mboost(data = test1$data,
                                 formula = TARGETVAR~bbs(WS100,knots=8),
                                 quantiles = seq(0.1,0.9,by=0.1),
                                 sort = T,
                                 sort_limits = list(U=0.999,L=0.001),
                                 cores = 8L)

# test with save model function 1 worker
tmp <- mqr_qreg_mboost(data = test1$data,
                       formula = TARGETVAR~bbs(WS100,knots=8),
                       quantiles = seq(0.1,0.9,by=0.1),
                       sort = T,
                       sort_limits = list(U=0.999,L=0.001),
                       save_models_path = "./tmp")
rm(tmp)

file.remove(list.files(pattern = "*.rda"))

# test with save model function 8 workers
tmp <- mqr_qreg_mboost(data = test1$data,
                       formula = TARGETVAR~bbs(WS100,knots=8),
                       quantiles = seq(0.1,0.9,by=0.1),
                       sort = T,
                       sort_limits = list(U=0.999,L=0.001),
                       save_models_path = "./tmp",
                       cores=8L)

file.remove(list.files(pattern = "*.rda"))
rm(tmp)





#####################

## back compatibility test
## test in series
test1$gam_mqr4 <- MQR_qreg_mboost(data = test1$data,
                                 formula = TARGETVAR~bbs(WS100,knots=8),
                                 bc_mstop = 100,
                                 bc_nu = 0.1,
                                 quantiles = seq(0.1,0.9,by=0.1),
                                 Sort = T,
                                 SortLimits = list(U=0.999,L=0.001))

# test with modified hyperparameters
test1$gam_mqr5 <- MQR_qreg_mboost(data = test1$data,
                                  formula = TARGETVAR~bbs(WS100,knots=8),
                                  bc_mstop = 100,
                                  bc_nu = 0.5,
                                  quantiles = seq(0.1,0.9,by=0.1),
                                  Sort = T,
                                  SortLimits = list(U=0.999,L=0.001))


# test in parallel ---> should equal model 1
test1$gam_mqr6 <- MQR_qreg_mboost(data = test1$data,
                                  formula = TARGETVAR~bbs(WS100,knots=8),
                                  bc_mstop = 100,
                                  bc_nu = 0.1,
                                  quantiles = seq(0.1,0.9,by=0.1),
                                  Sort = T,
                                  SortLimits = list(U=0.999,L=0.001),
                                  cores = 8L,
                                  parallel = TRUE)



### default
test1$gam_mqr7 <- mqr_qreg_mboost(data = test1$data,
                                  formula = TARGETVAR~WS100,
                                  quantiles = seq(0.1,0.9,by=0.1),
                                  sort = T,
                                  sort_limits = list(U=0.999,L=0.001),
                                  cores = 8L,
                                  baselearner="btree")








par(mar=c(3,3,0.5,1))  # Trim margin around plot [b,l,t,r]
par(tcl=0.35)  # Switch tick marks to insides of axes
par(mgp=c(1.5,0.2,0))  # Set margin lines; default c(3,1,0) [title,labels,line]
par(xaxs="r",yaxs="r")  # Extend axis limits by 4% ("i" does no extension)

i_ts <- unique(test1$data$ISSUEdtm)[3]

plot(test1$gam_mqr[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

plot(test1$gam_mqr3[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

plot(test1$gam_mqr4[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

plot(test1$gam_mqr6[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)

plot(test1$gam_mqr7[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)



plot(test1$gam_mqr2[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)


plot(test1$gam_mqr5[which(test1$data$ISSUEdtm==i_ts),],xlab="Time Index [Hours]",ylab="Power [Capacity Factor]",axes=F,Legend = 1,ylim=c(0,1)); axis(1,1:24,pos=-0.07); axis(2,las=1)
lines(test1$data$TARGETVAR[which(test1$data$ISSUEdtm==i_ts)],lwd=3)



lapply(names(test1)[-1][c(1,3,4,6,7,2,5)],function(x){
  
  reliability(qrdata = test1[[x]],
              realisations = test1$data$TARGETVAR,
              kfolds = test1$data$kfold,main=x)
  
  
})


lapply(names(test1)[-1][c(1,3,4,6,7,2,5)],function(x){
  
  pinball(qrdata = test1[[x]],
          realisations = test1$data$TARGETVAR,
          kfolds = test1$data$kfold,main=x)
  
  
})




#################################################################################
### Parametric Regression using GAMLSS ####
#################################################################################



load_all("./", TRUE)

test1<-list(data=Wind)

# library(gamlss)
## test in series
test1$ppd <- ppd_gamlss(data = test1$data,
                         formula = TARGETVAR~bs(WS100,df=3),
                         sigma.formula = ~WS100,
                         nu.formula = ~WS100,
                         tau.formula = ~WS100,
                         family =  gamlss.dist::BEINF(),
                         method = mixed(20,10))

summary(test1$ppd$Test)



# test with modified hyperppdmeters
test1$ppd2 <- ppd_gamlss(data = test1$data,
                         formula = TARGETVAR~bs(WS100,df=3),
                         sigma.formula = ~WS100,
                         nu.formula = ~WS100,
                         tau.formula = ~WS100,
                         family =  gamlss.dist::BEINF(),
                         sigma.start = 0.05,#NO,  #
                         method = mixed(20,10))

summary(test1$ppd2$Test)

identical(coef(test1$ppd2$Test),coef(test1$ppd$Test))



# test in pllel ---> should equal model 1
test1$ppd3 <- ppd_gamlss(data = test1$data,
                         formula = TARGETVAR~bs(WS100,df=3),
                         sigma.formula = ~WS100,
                         nu.formula = ~WS100,
                         tau.formula = ~WS100,
                         family =  gamlss.dist::BEINF(),
                         method = mixed(20,10),
                         cores = 6)

summary(test1$ppd3$Test)

identical(coef(test1$ppd3$Test),coef(test1$ppd$Test))



# test with save model function 1 worker
tmp <- ppd_gamlss(data = test1$data,
                  formula = TARGETVAR~bs(WS100,df=3),
                  sigma.formula = ~WS100,
                  nu.formula = ~WS100,
                  tau.formula = ~WS100,
                  family =  gamlss.dist::BEINF(),
                  method = mixed(20,10),
                  save_models_path = "./tmp")
rm(tmp)

file.remove(list.files(pattern = "*.rda"))

# test with save model function 8 workers
tmp <-  ppd_gamlss(data = test1$data,
                   formula = TARGETVAR~bs(WS100,df=3),
                   sigma.formula = ~WS100,
                   nu.formula = ~WS100,
                   tau.formula = ~WS100,
                   family =  gamlss.dist::BEINF(),
                   method = mixed(20,10),
                   cores=8L,
                   save_models_path = "./tmp")

file.remove(list.files(pattern = "*.rda"))
rm(tmp)





#####################

## back compatibility test
## test in series
# library(gamlss)
test1$ppd4 <- Para_gamlss(data = test1$data,
                          formula = TARGETVAR~bs(WS100,df=3),
                          sigma.formula = ~WS100,
                          nu.formula = ~WS100,
                          tau.formula = ~WS100,
                          family =  gamlss.dist::BEINF(),
                          method = mixed(20,10))

identical(coef(test1$ppd4$Test),coef(test1$ppd$Test))



# test with modified hyperparameters
test1$ppd5 <- Para_gamlss(data = test1$data,
                          formula = TARGETVAR~bs(WS100,df=3),
                          sigma.formula = ~WS100,
                          nu.formula = ~WS100,
                          tau.formula = ~WS100,
                          family =  gamlss.dist::BEINF(),
                          sigma.start=0.05,
                          method = mixed(20,10))

identical(coef(test1$ppd5$Test),coef(test1$ppd2$Test))


# test in parallel ---> should equal model 1
test1$ppd6 <- Para_gamlss(data = test1$data,
                          formula = TARGETVAR~bs(WS100,df=3),
                          sigma.formula = ~WS100,
                          nu.formula = ~WS100,
                          tau.formula = ~WS100,
                          family =  gamlss.dist::BEINF(),
                          method = mixed(20,10),
                          parallel = TRUE,
                          cores = 6)

identical(coef(test1$ppd6$Test),coef(test1$ppd$Test))


### default


test1$ppd7 <- ppd_gamlss(data = test1$data,
                         formula = TARGETVAR~bs(WS100,df=3),
                         sigma.formula = ~WS100,
                         nu.formula = ~WS100,
                         tau.formula = ~WS100,
                         family =  gamlss.dist::BEINF())

test1$ppd8 <- Para_gamlss(data = test1$data,
                          formula = TARGETVAR~bs(WS100,df=3),
                          sigma.formula = ~WS100,
                          nu.formula = ~WS100,
                          tau.formula = ~WS100,
                          family =  gamlss.dist::BEINF())



identical(coef(test1$ppd7$Test),coef(test1$ppd8$Test))





#################################################################################
### Parametric Regression using gamboostLSS ####
#################################################################################





load_all("./", TRUE)

test1<-list(data=Wind)


### when the function runs in series, each kfold model is started from the same offset in default mode ---> family = gamboostLSS::GaussianLSS().
### this looks as if it's due to the "family" argument retaining a default after it's ran once
### in parallel the mean of each subset is used as the offset, which is what we would expect.
### this isn't a big deal really, except for **exact** reproducability when changing from parallel to series computation
### let's just override the default behavior and start all models at the same offset in case we want to turn this into a real test

test1$ppd <- ppd_gamboostlss(data = test1$data,
                             formula = TARGETVAR~bbs(WS100),
                             family = gamboostLSS::GaussianLSS(mu=.3,sigma = .3))

summary(test1$ppd$Test)


# test with modified hyperppdmeters
test1$ppd2 <- ppd_gamboostlss(data = test1$data,
                              formula = TARGETVAR~bbs(WS100),
                              control = mboost::boost_control(mstop=100,nu=0.25),
                              family = gamboostLSS::GaussianLSS(mu=.3,sigma = .3))



identical(coef(test1$ppd2$Test),coef(test1$ppd$Test))


# test in pllel ---> should equal model 1 --- 
test1$ppd3 <- ppd_gamboostlss(data = test1$data,
                              formula = TARGETVAR~bbs(WS100),
                              cores = 4,
                              family = gamboostLSS::GaussianLSS(mu=.3,sigma = .3))


identical(coef(test1$ppd3$Test),coef(test1$ppd$Test))


# test with save model function 1 worker
tmp <- ppd_gamboostlss(data = test1$data,
                       formula = TARGETVAR~bbs(WS100),
                       family = gamboostLSS::GaussianLSS(mu=.3,sigma = .3),
                       save_models_path = "./tmp")
rm(tmp)

file.remove(list.files(pattern = "*.rda"))

# test with save model function 8 workers
tmp <-  ppd_gamboostlss(data = test1$data,
                        formula = TARGETVAR~bbs(WS100),
                        family = gamboostLSS::GaussianLSS(mu=.3,sigma = .3),
                        save_models_path = "./tmp",
                        cores=8L)

file.remove(list.files(pattern = "*.rda"))
rm(tmp)





#####################

## back compatibility test
## test in series
# library(gamlss)
test1$ppd4 <- Para_gamboostLSS(data = test1$data,
                               formula = TARGETVAR~bbs(WS100),
                               families = gamboostLSS::GaussianLSS(mu=.3,sigma = .3),
                               parallel = FALSE)

identical(coef(test1$ppd4$Test),coef(test1$ppd$Test))



# test with modified hyperparameters
test1$ppd5 <- Para_gamboostLSS(data = test1$data,
                               formula = TARGETVAR~bbs(WS100),
                               families = gamboostLSS::GaussianLSS(mu=.3,sigma = .3),
                               parallel = FALSE,
                               control = mboost::boost_control(mstop=100,nu=0.25))

identical(coef(test1$ppd5$Test),coef(test1$ppd2$Test))


# test in parallel ---> should equal model 1
test1$ppd6 <- Para_gamboostLSS(data = test1$data,
                               formula = TARGETVAR~bbs(WS100),
                               families = gamboostLSS::GaussianLSS(mu=.3,sigma = .3),
                               parallel = TRUE,
                               cores = 6)




identical(coef(test1$ppd6$Test),coef(test1$ppd$Test))


### default

test1$ppd7 <- ppd_gamboostlss(data = test1$data,
                              formula = TARGETVAR~bbs(WS100))

test1$ppd8 <- Para_gamboostLSS(data = test1$data,
                               formula = TARGETVAR~bbs(WS100))



identical(coef(test1$ppd7$Test),coef(test1$ppd8$Test))








########################################################################################

# example of offset problem --- run 1st half of fitting function then this

# 
# for (i in unique(tempdata$kfold)){
#   
#   
#   
#   GAMLSSmodelList[[i]] <- gamboostLSS::gamboostLSS(data = tempdata[tempdata$kfold!=i & tempdata$kfold != "Test",],
#                                    formula = formula,
#                                    families = family)
#   
#   
#   
#   
#   
# }
# 
# ### this should be the mean of the subset of data as below! But they're all the same
# GAMLSSmodelList$`Fold 1`$mu$offset
# GAMLSSmodelList$Test$mu$offset
# 
# mean(tempdata[tempdata$kfold!="Fold 1" & tempdata$kfold != "Test",]$TARGETVAR,na.rm=T)
# mean(tempdata[tempdata$kfold != "Test",]$TARGETVAR,na.rm=T)
#
# the models contain all the correct subsets of data though, think it comes down to the <<family>> object we defined updating with a default...








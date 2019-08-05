### StrathCast Extended Example
require(devtools)
require(roxygen2)
require(rstudioapi)

PackagePath <- dirname(getActiveDocumentContext()$path)
setwd(PackagePath)

# Install from local repository
install(PackagePath)
# Update package documentation
document(pkg = ".")
# Load Package
require(ProbCast)



setwd(PackagePath)
### Testing Functionalisty of StrathCast #####

## Add some features first...

Wind$WS100 <- sqrt(Wind$U100^2+Wind$V100^2)
Wind$Power <- pmin(Wind$WS100,11)^3

## Set-up simple kfold CV
Wind$kfold <- c(rep(c("fold1","fold2"),each=6000),rep("Test",nrow(Wind)-12000))
# Wind$kfold <- c(rep(c("fold1","fold2","fold3","fold4"),each=nrow(Wind)/4))

### Multiple Quantile Regression using GBM ####
test1<-list(data=Wind)

test1$gbm_mqr <- MQR_gbm(data = test1$data,
                         formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                         gbm_params = list(interaction.depth = 3,
                                                        n.trees = 1000,
                                                        shrinkage = 0.05,
                                                        cv.folds = 0,
                                                        n.minobsinnode = 20,
                                                        bag.fraction = 0.5,
                                                        keep.data = F),
                         quantiles = seq(0.1,0.9,by=0.1),
                         Sort = T,
                         SortLimits = list(U=0.999,L=0.001))

test1$gbm_mqr <- MQR_gbm(data = test1$data,
                         formula = TARGETVAR~U100+V100+U10+V10+(sqrt((U100^2+V100^2))),
                         gbm_params = list(interaction.depth = 3,
                                           n.trees = 1000,
                                           shrinkage = 0.05,
                                           cv.folds = 0,
                                           n.minobsinnode = 20,
                                           bag.fraction = 0.5,
                                           keep.data = F),
                         parallel = T,
                         cores = 3,
                         quantiles = seq(0.1,0.9,by=0.1),
                         Sort = T,
                         SortLimits = list(U=0.999,L=0.001))


plot(test1$gbm_mqr[1:240,],xlab="Time Index",ylab="Power")
lines(test1$data$TARGETVAR[1:240])

reliability(qrdata = test1$gbm_mqr,
            realisations = test1$data$TARGETVAR,
            kfolds = test1$data$kfold)

pinball(qrdata = test1$gbm_mqr,
         realisations = test1$data$TARGETVAR,
         kfolds = test1$data$kfold)

index <- 1000
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],method = "spline")
plot(cdf(seq(0,1,by=0.001)),type="l")
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],method = "linear")
lines(cdf(seq(0,1,by=0.001)),lty=2,col=2)
cdf <- contCDF(quantiles = test1$gbm_mqr[index,],kfold = NA,method = "spline", tails=list(method="exponential",L=0,U=1,nBins=5,preds=test1$gbm_mqr,targetvar=test1$data$TARGETVAR,ntailpoints=25))
lines(cdf(seq(0,1,by=0.001)),lty=3,col=3)


# test1$X_gbm <- PIT(test1$gbm_mqr,test1$data$TARGETVAR,method = "linear",tails=list(method="exponential",L=0,U=1,nBins=5,preds=test1$gbm_mqr,targetvar=test1$data$TARGETVAR,ntailpoints=25))
test1$X_gbm <- PIT(test1$gbm_mqr,test1$data$TARGETVAR,method = "linear",tails=list(method="interpolate",L=0,U=1))
hist(test1$X_gbm,breaks = 100,freq=F)

### Parametric PredDist Using GAMLSS ####
# require(gamlss.tr)
# gen.trun(family = NO,par = c(0,1),type = "both",name="test")


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

test1$gamlss_mqr <- PPD_2_MultiQR(data=test1$data,
                                  models = test1$ppd,
                                  params = F)


reliability(qrdata = test1$gamlss_mqr,
            realisations = test1$data$TARGETVAR,
            kfolds = test1$data$kfold)

pinball(qrdata = test1$gamlss_mqr,
         realisations = test1$data$TARGETVAR,
         kfolds = test1$data$kfold)


# test1$data[test1$data$TARGETVAR<0 | test1$data$TARGETVAR>1,]
test1$X_gamlss <- PIT(test1$ppd,data = test1$data)
hist(test1$X_gamlss,breaks = 100,freq = F,ylim = c(0,2)); lines(c(0,1),c(1,1),lty=2)



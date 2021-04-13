### StrathCast Extended Example
require(devtools)
require(roxygen2)
require(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))


PackagePath <- ".."
# setwd(PackagePath)

# Update package documentation
document(pkg = PackagePath)
# Install from local repository
install(PackagePath)
# Load Package
require(ProbCast)


### Testing Functionalisty of StrathCast #####

## Add some features first...
Wind <- data.table(Wind)
Wind$WS100 <- sqrt(Wind$U100^2+Wind$V100^2)
Wind$WS10 <- sqrt(Wind$U10^2+Wind$V10^2)
Wind$Power <- pmin(Wind$WS100,11)^3

## Set-up simple kfold CV. NB --- For scenario forecasting make sure the CV folds don't cross issue times
Wind$kfold <- "fold 1"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2012-06-30",tz="UTC")] <- "fold 2"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2012-12-31",tz="UTC")] <- "fold 3"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2013-06-30",tz="UTC")] <- "Test"
Wind$BadData <- runif(n = nrow(Wind))<0.005

QR_form <- ~ gam_pred + WS100 + WS10
R2_form <- ~ WS100 + WS10

## GAM model with some quantiles
Model_1 <- qreg_gam(data = Wind,
                    formula = TARGETVAR ~ te(WS100,WS10,k=5) + te(U100,V100,k=12),
                    formula_qr = NULL,
                    # formula_qr = ~ gam_pred + WS100 + WS10,
                    # formula_qr =  QR_form,
                    # cv_folds = NULL,
                    cv_folds = "kfold",
                    # cv_folds = Wind$kfold,
                    # model_res2 = T,
                    # formula_res2 = R2_form,
                    # exclude_train = NULL,
                    # exclude_train = "BadData",
                    exclude_train = Wind$BadData,
                    quantiles = c(0.05,0.25,0.5,0.75,0.95),sort_limits = list(U=1,L=0),
                    discrete=F)


summary(Model_1)


## Check GAM fits...
# gam.check(Model_1$models$gams$`Fold 1`)
# gam.check(Model_1$models$gams$`Fold 1_r`)


## Plot...
{indexes <- 1:240 + 240*4
  plot(Model_1$mqr_pred[indexes,],targetTimes = Wind$TARGETdtm[indexes])
  lines(Wind$TARGETdtm[indexes],Wind$TARGETVAR[indexes],lwd=2)}


## Add some extra quantiles
Model_1 <- qreg_gam.add_quantiles(Model_1,data = Wind,quantiles = c(0.1,0.9))

{indexes <- 1:240 + 240*4
  plot(Model_1$mqr_pred[indexes,],targetTimes = Wind$TARGETdtm[indexes])
  lines(Wind$TARGETdtm[indexes],Wind$TARGETVAR[indexes],lwd=2)}

## Test predict S3 method
Temp_predict <- predict(Model_1,newdata=Wind[kfold=="Test",],
                        sort_limits = list(U=1,L=0))

plot(Temp_predict$mqr_pred[1:240,],q50_line = T)
reliability(qrdata = Model_1$mqr_pred,realisations = Wind$TARGETVAR,kfolds = Wind$kfold)


## Test update function
Model_1_ud <- qreg_gam.update(object = Model_1,newdata = Wind[ISSUEdtm>as.POSIXct("2013-06-30",tz="UTC") & 
                                                                ISSUEdtm<as.POSIXct("2013-08-01",tz="UTC"),])

Temp_predict <- predict(Model_1_ud,newdata=Wind[ISSUEdtm>=as.POSIXct("2013-08-01",tz="UTC"),],
                        sort_limits = list(U=1,L=0))

plot(Temp_predict$mqr_pred[1:240,],q50_line = T)



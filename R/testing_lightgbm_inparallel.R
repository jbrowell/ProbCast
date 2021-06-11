require(data.table)
#require(graphics)
require(ProbCast)
require(gbm)
require(rstudioapi)

# setwd(dirname(getActiveDocumentContext()$path))
setwd("C:/Users/Ciaran/Desktop/")
## testing out putting lightgbm in parallel.. what is making it fail?


#setwd(dirname(getActiveDocumentContext()$path))
cores <- detectCores()-2
md <- 5
nl <- 15
nfolds <- 5
formula = TARGETVAR~U100+V100+U10+V10+WS100
features <- labels(terms(formula))
response <- as.character(formula)[[2]]


quantiles = seq(0.05,0.95,by=0.05)

Wind$WS100 <- sqrt(Wind$U100^2+Wind$V100^2)
Wind$Power <- pmin(Wind$WS100,11)^3

Wind$kfold <- "Fold 1"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2012-06-30",tz="UTC")] <- "Fold 2"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2012-12-31",tz="UTC")] <- "Fold 3"
Wind$kfold[Wind$ISSUEdtm>as.POSIXct("2013-06-30",tz="UTC")] <- "Test"
#Wind <- as.data.table(Wind)


# set up parallel workers
cl <- parallel::makeCluster(cores)
doSNOW::registerDoSNOW(cl)
#set up progress bar
iterations <- length(quantiles)
pb <- utils::txtProgressBar(max = iterations, style = 3)
progress <- function(n) utils::setTxtProgressBar(pb, n)
opts <- list(progress = progress)
gc()

output <- foreach::foreach(q = quantiles,.packages ="lightgbm",.options.snow = opts) %dopar% {
  models <- list()
  for (fold in unique(Wind$kfold)){
    train_data <- Wind[Wind$kfold != fold & Wind$kfold != "Test" & !is.na(Wind[[formula[[2]]]]),]
    X <- as.matrix(train_data[, features])
    dataset <- lgb.Dataset(data=X, label=train_data[,response])
    #lgb_params <- list(..., num_iterations=num_iterations, objective="quantile")
    lgb_params <- list(objective="quantile", num_iterations=100)
    # lgb_params <- list(objective="regression", num_iterations=100)
    lgb_model = lgb.train(params=lgb_params, data=dataset, alpha=q, num_threads=1)
    # lgb_model = lgb.train(params=lgb_params, data=dataset, num_threads=1)
    model_file = paste0(fold,"_",q,".txt")
    lgb.save(lgb_model, model_file)
    models[[fold]] <- lgb_model
    
    #test_data <- Wind[Wind$kfold==fold, features]
    #models[[paste0(fold,"_fc")]] <- predict(lgb_model, as.matrix(test_data))
    
    
  }
  
  return(models)
  
}

close(pb)
parallel::stopCluster(cl)
names(output) <- paste0(quantiles)

#  crashing just like before!
# test_data <- as.matrix(Wind[Wind$kfold=="Fold 3", features])
# lightgbm:::predict.lgb.Booster(output$`0.05`$`Fold 3`, test_data,num_iteration = 100)


load_booster <- lightgbm::lgb.load(filename = "Fold 3_0.05.txt")
test_data <- as.matrix(Wind[Wind$kfold=="Fold 3", features])
plot(head(lightgbm:::predict.lgb.Booster(load_booster, test_data),48),ylim=c(0,1),type="l")


load_booster <- lightgbm::lgb.load(filename = "Fold 3_0.95.txt")
test_data <- as.matrix(Wind[Wind$kfold=="Fold 3", features])
lines(lightgbm:::predict.lgb.Booster(load_booster, test_data),col="red")




########################################## this works for me
cl <- parallel::makeCluster(cores)
doSNOW::registerDoSNOW(cl)
#set up progress bar
iterations <- length(quantiles)
pb <- utils::txtProgressBar(max = iterations, style = 3)
progress <- function(n) utils::setTxtProgressBar(pb, n)
opts <- list(progress = progress)
gc()

forecasts <- foreach::foreach(q = quantiles,.packages ="lightgbm",.options.snow = opts) %dopar% {
  foldresults <- list()
  ## now make predictions..
  for (fold in unique(Wind$kfold)){
    # test_data <- Wind[Wind$kfold==fold, features]
    # foldresults[[fold]] <- try(lightgbm:::predict.lgb.Booster(output[[q]][[fold]], as.matrix(test_data)))
    
    
    load_booster <- lgb.load(filename = paste0(fold,"_",q,".txt"))
    test_data <- as.matrix(Wind[Wind$kfold==fold, features])
    foldresults[[fold]] <- predict(load_booster, test_data)
    
    
    
  }
  
  return(foldresults)
}
close(pb)
parallel::stopCluster(cl)
names(forecasts) <- paste0(quantiles)



output  <- data.frame(matrix(NA,ncol = length(quantiles), nrow = nrow(Wind)))
colnames(output) <- paste0("q",100*quantiles)

for(q in quantiles){
  for(fold in unique(Wind$kfold)){
    
    output[Wind$kfold==fold,paste0("q",q*100)] <-forecasts[[as.character(q)]][[fold]]
    
  }
}
# class of new cv object
class(output) <- c("MultiQR","data.frame")


plot(output[1:100,])





# here's my session info --->
sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19041)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United Kingdom.1252 
# [2] LC_CTYPE=English_United Kingdom.1252   
# [3] LC_MONETARY=English_United Kingdom.1252
# [4] LC_NUMERIC=C                           
# [5] LC_TIME=English_United Kingdom.1252    
# 
# attached base packages:
#   [1] parallel  splines   stats     graphics  grDevices utils    
# [7] datasets  methods   base     
# 
# other attached packages:
#   [1] rstudioapi_0.11     gbm_2.1.8           ProbCast_0.0.0.9000
# [4] quantreg_5.75       SparseM_1.78        mgcv_1.8-31        
# [7] evgam_0.1.4         fst_0.9.4           mvnfast_0.2.5.1    
# [10] gamboostLSS_2.0-1.1 mboost_2.9-3        stabs_0.6-3        
# [13] doSNOW_1.0.19       snow_0.4-3          iterators_1.0.13   
# [16] foreach_1.5.1       gamlss_5.2-0        nlme_3.1-148       
# [19] gamlss.dist_5.1-7   MASS_7.3-51.6       gamlss.data_5.1-4  
# [22] data.table_1.14.0  
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.6         compiler_4.0.2     tools_4.0.2       
# [4] partykit_1.2-10    rpart_4.1-15       jsonlite_1.7.0    
# [7] lattice_0.20-41    Matrix_1.2-18      mvtnorm_1.1-1     
# [10] libcoin_1.0-6      MatrixModels_0.4-1 grid_4.0.2        
# [13] R6_2.4.1           survival_3.1-12    lightgbm_3.2.1    
# [16] Formula_1.2-4      conquer_1.0.2      matrixStats_0.57.0
# [19] codetools_0.2-16   nnls_1.4           quadprog_1.5-8    
# [22] inum_1.0-1    










fold <- unique(Wind$kfold)[1]
test_data <- Wind[Wind$kfold==fold, features]
fcs <- predict(output[[0.05]][[fold]], as.matrix(test_data))

## get error "no applicable method for 'predict' applied to an object of class "NULL"
fold <- unique(Wind$kfold)[1]
q <- quantiles[1]
class(output[[fold]][[q]])

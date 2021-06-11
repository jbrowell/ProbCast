require(data.table)
#require(graphics)
require(ProbCast)
require(gbm)
require(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))
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
    #lgb_params <- list(objective="quantile", num_iterations=100)
    lgb_params <- list(objective="regression", num_iterations=100)
    #lgb_model = lgb.train(params=lgb_params, data=dataset, alpha=q, num_threads=1)
    lgb_model = lgb.train(params=lgb_params, data=dataset, num_threads=1)
    model_file = tempfile(fileext=paste0(fold, ".txt"))
    lgb.save(lgb_model, model_file)
    #models[[fold]] <- lgb_model
    
    #test_data <- Wind[Wind$kfold==fold, features]
    #models[[paste0(fold,"_fc")]] <- predict(lgb_model, as.matrix(test_data))
    
    
  }
  
  return (models)
  
}

close(pb)
parallel::stopCluster(cl)
names(output) <- paste0(quantiles)


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
    test_data <- Wind[Wind$kfold==fold, features]
    foldresults[[fold]] <- try(lightgbm:::predict.lgb.Booster(output[[q]][[fold]], as.matrix(test_data)))
  }
  
  #test <- unlist(foldresults)
  #return(test)
}
close(pb)
parallel::stopCluster(cl)

fold <- unique(Wind$kfold)[1]
test_data <- Wind[Wind$kfold==fold, features]
fcs <- predict(output[[0.05]][[fold]], as.matrix(test_data))

## get error "no applicable method for 'predict' applied to an object of class "NULL"
fold <- unique(Wind$kfold)[1]
q <- quantiles[1]
class(output[[fold]][[q]])

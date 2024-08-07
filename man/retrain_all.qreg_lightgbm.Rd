% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MQR_lightgbm.R
\name{retrain_all.qreg_lightgbm}
\alias{retrain_all.qreg_lightgbm}
\title{Method for repeated model retraining throughout the Test set, based on  Multiple 
Quantile Regression from lightgbm Gradient Boosted Decision Trees.}
\usage{
\method{retrain_all}{qreg_lightgbm}(
  object,
  data,
  retrain_daily_frequency,
  issue_datetime_column,
  cv_folds,
  pckgs = NULL,
  cores = 1,
  ...
)
}
\arguments{
\item{object}{fitted model of class \code{qreg_lightgbm} obtained from the function \code{qreg_lightgbm()}.}

\item{data}{data.frame of model inputs and observations, including a column containing fold labels: 
models will be retrained throughout the subset of the data labelled 'Test'.}

\item{retrain_daily_frequency}{Frequency, in number of days, which the model will be retrained over the Test set.}

\item{issue_datetime_column}{Name of the column in data containing the forecast issue times. Column of issue times must be POSIXct or Date.}

\item{cv_folds}{Name of column in data containing the fold labels. Models will 
only be retrained for time points beyond the start of the Test fold.}

\item{pckgs}{specify additional packages required for
each worker (e.g. c("data.table") if data stored as such).}

\item{cores}{the number of available cores. Defaults to one, i.e. no parallelisation, although in this case the user
must still specify \code{pckgs} if applicable.}

\item{...}{additional arguments; not currently used.}
}
\value{
a named list containing fitted models as a list of \code{qreg_lightgbm} objects.
}
\description{
This function returns updated fitted \code{qreg_lightgbm} models and predictions, which are the 
supplied previously fitted models updated throughout the Test set. S3 Method for \code{qreg_lightgbm} objects
}
\details{
this function periodically retrains the fitted lightgbm model over the Test set and updates the predictions. 
NOTE predictions saved within the lightgbm object are not updated within this function.
}
\author{
Gordon McFadzean, \email{gordon.mcfadzean@tneigroup.com}; Rosemary Tawn, \email{rosemary.tawn@tneigroup.com}
}
\keyword{Quantile}
\keyword{Regression}

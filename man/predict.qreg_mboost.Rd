% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MQR_qreg_mboost.R
\name{predict.qreg_mboost}
\alias{predict.qreg_mboost}
\title{Predict method for Multiple Quantile Regression using \code{mboost}}
\usage{
\method{predict}{qreg_mboost}(
  object,
  newdata = NULL,
  quantiles = NULL,
  model_name = NULL,
  sort = T,
  sort_limits = NULL,
  ...
)
}
\arguments{
\item{object}{object of class \code{qreg_mboost} obtained from the function \code{qreg_mboost()}.}

\item{newdata}{data.frame of observations for which to make predictions}

\item{quantiles}{the quantiles to predict. Default is all the quantiles present in \code{object}}

\item{model_name}{the name of the model to be used for prediction.
Unless specified, \code{object$default_model} is used.}

\item{sort}{sort quantiles using \code{SortQuantiles()}?}

\item{sort_limits}{\code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and 
lower limits given by \code{list(U=upperlim,L=lowerlim)}.}

\item{...}{additional arguments passed on \code{mboost::predict.mboost()}}
}
\value{
Quantile forecasts in a \code{MultiQR} object.
}
\description{
This function returns multiple quantile predictions as an object of class \code{MultiQR}
based on a ProbCast mboost fitted model. S3 Method for for \code{qreg_mboost} objects
}
\details{
this function returns predictive quantiles for each row in \code{newdata},
the result is returned as a \code{MultiQR} object
}
\author{
Ciaran Gilbert, \email{ciaran.gilbert@strath.ac.uk}
}
\keyword{Quantile}
\keyword{Regression}

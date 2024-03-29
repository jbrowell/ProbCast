% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MQR_qreg_gam.R
\name{predict.qreg_gam}
\alias{predict.qreg_gam}
\title{Predict from model based on Generalised Additive Model and Linear Quantile Regression}
\usage{
\method{predict}{qreg_gam}(
  object,
  newdata = NULL,
  quantiles = NULL,
  model_name = NULL,
  sort = T,
  sort_limits = NULL
)
}
\arguments{
\item{object}{An \code{qreg_gam} object containing the model to predict from.}

\item{newdata}{A data frame or data table containing the values of the model
covariates at which predictions are required.}

\item{quantiles}{The probability levels at which quantile predictions should
be produced.}

\item{model_name}{The name of the model in \code{object} to be used for prediction. E.g.
Specific cross-vlaidation fold, test model or some other version.}

\item{sort}{\code{boolean} Sort quantiles using \code{SortQuantiles()}?}

\item{sort_limits}{\code{Limits} argument to be passed to \code{SortQuantiles()}. 
Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.}
}
\value{
A list with elements \code{gam_pred}, deterministic predictions (conditional expectation) from
main GAM model, and \code{mqr_pred}, multiple predictive quantiles in a \code{MultiQR} object.
}
\description{
This function predicts from multiple conditional linear quantile regression models of the residuals of
a generalised additive model fit using \code{qreg_gam}.
}
\details{
Predict method for multiple quantile regression models of the class \code{qreg_gam}.
}
\author{
Jethro Browell, \email{jethro.browell@glasgow.ac.uk}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MQR_rq.R
\name{qreg_mrq}
\alias{qreg_mrq}
\title{Multiple quantile regression, and quantile dressing, using linear quantile regression}
\usage{
qreg_mrq(
  data,
  formula,
  quantiles = c(0.25, 0.5, 0.75),
  offset = NULL,
  cv_folds = NULL,
  exclude_train = NULL,
  sort = T,
  sort_limits = NULL,
  ...
)
}
\arguments{
\item{data}{A \code{data.frame} containing target and explanatory
variables. Optionally, supplying a deterministic forecast with \code{offset} will return
the this forecast dressed with multiple predictive quantiles.
May optionally contain a column called "kfold" with numbered/labelled folds and "Test" for test data.}

\item{quantiles}{The quantiles to fit models for.}

\item{offset}{The column name in \code{data} of an optional deterministic forecast to be dressed with multiple predictive quantiles.}

\item{cv_folds}{Control for cross-validation with various options, either:
\itemize{
 \item the column name of the fold index supplied in data. Observations and inputs 
  in the index labelled "Test" will serve as test data and held out in model training.
 \item an integer giving the number of cross validation folds to generate. Folds are constructed as block chunks. 
 Default behaviour is 5 folds.
 \item NULL indicates that no cross validation should be performed and the returned model is trained on all \code{data}.
}}

\item{exclude_train}{A column name in \code{data} indicating if a row should be excluded from model
training, i.e. if it contains bad data (will be coerced to \code{logical}). Alternatively,
an \code{integer} or \code{logical} vector with length equal to the number of rows in \code{data} indicating
the same. Rows labelled \code{TRUE} are excluded from model training.}

\item{sort}{\code{boolean} Sort quantiles using \code{SortQuantiles()}?}

\item{sort_limits}{\code{Limits} argument to be passed to \code{SortQuantiles()}. Constrains quantiles to upper and lower limits given by \code{list(U=upperlim,L=lowerlim)}.}

\item{...}{Additional arguments passed to \code{rq}.}

\item{formala}{A \code{formula} object with the response on the left
of an ~ operator, and the terms, separated by + operators, on the right passed to \code{rq()}.}
}
\value{
Returns a \code{list} containing predictive quantiles (in a \code{MultiQR} object) and \code{rq} models.
}
\description{
This function fits multiple conditional linear quantile regression models, optionally to the residuals of
a user-specified deterministic forecast with facilities for cross-validation.
}
\details{
The returned predictive quantiles are those produced out-of-sample for each
cross-validation fold (using models trained on the remaining folds but not "Test" data).
Predictive quantiles corresponding to "Test" data are produced using models trained on all
non-test data.
}
\author{
Jethro Browell, \email{jethro.browell@glasgow.ac.uk}
}
\keyword{Quantile}
\keyword{Regression}

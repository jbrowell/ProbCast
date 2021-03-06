% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{cv_control}
\alias{cv_control}
\title{cross validation set-up for model fitting functions}
\usage{
cv_control(data, cv_folds = "kfold")
}
\arguments{
\item{data}{A \code{data.frame} containing target and explanatory
variables. May optionally contain a column with labeled folds and "Test" for test data. See \code{cv_folds}.}

\item{cv_folds}{Control for cross-validation with various options, either:
\itemize{
 \item the column name of the fold index supplied in data. Using this option, if there is a 
  fold called "Test", this will serve as test data and held out in the model training.
 \item an integer giving the number of cv evaluations to perform. Folds are constructed as block chunks. 
 Default behaviour is 5 folds.
 \item NULL indicates that no cross validation should be performed and the returned model is trained on all \code{data}.
 \item vector of length==nrow(data) containing character or numeric fold labels
}}
}
\value{
A named list containing the kfold indexes and the kfold loop reference values
}
\description{
This generates kfold indexes and vectors of reference values for the model fitting loop
}
\details{
....
}
\author{
Ciaran Gilbert, \email{ciaran.gilbert@strath.ac.uk}
}

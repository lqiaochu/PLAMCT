% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plamct.R
\name{resample}
\alias{resample}
\title{Title resample}
\usage{
resample(weight, res, fit)
}
\arguments{
\item{weight}{weight}

\item{res}{residuals}

\item{fit}{fitted value}
}
\value{
vector ynew
}
\description{
Title resample
}
\examples{
y.boot=apply(w.boot,2,resample,res=residuals.b,fit=result$fitted)
}

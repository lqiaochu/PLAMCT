% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plamct.R
\name{regre}
\alias{regre}
\title{Title regre}
\usage{
regre(y, design_x)
}
\arguments{
\item{y}{response}

\item{design_x}{design matrix}
}
\value{
list
}
\description{
Title regre
}
\examples{
result.boot = apply(y.boot,2,regre,design_x=design_x)
}

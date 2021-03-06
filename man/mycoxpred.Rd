% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survfit.coxnet.R
\name{mycoxpred}
\alias{mycoxpred}
\title{Helper function to amend ... for new data in survfit.coxnet}
\usage{
mycoxpred(object, s, ...)
}
\arguments{
\item{object}{A class \code{coxnet} object.}

\item{s}{The response for the fitted model.}

\item{...}{The same ... that was passed to survfit.coxnet.}
}
\description{
This function amends the function arguments passed to survfit.coxnet
via ... if new data was passed to survfit.coxnet. It's a separate
function as we have to deal with function options like newstrata
and newoffset.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survfit.cv.glmnet.R
\name{survfit.cv.glmnet}
\alias{survfit.cv.glmnet}
\title{Compute a survival curve from a cv.glmnet object}
\usage{
\method{survfit}{cv.glmnet}(formula, s = c("lambda.1se", "lambda.min"), ...)
}
\arguments{
\item{formula}{A class \code{cv.glmnet} object. The object should have
been fit with \code{family = "cox"}.}

\item{s}{Value(s) of the penalty parameter lambda at which predictions
are required. Default is the value s="lambda.1se" stored on the CV object.
Alternatively s="lambda.min" can be used. If s is numeric, it is taken
as the value(s) of lambda to be used.}

\item{...}{Other arguments to be passed to \code{survfit.coxnet}.}
}
\value{
If \code{s} is a single value, an object of class "survfitcox"
and "survfit" containing one or more survival curves. Otherwise, a list
of such objects, one element for each value in \code{s}.
Methods defined for survfit objects are print, summary and plot.
}
\description{
Computes the predicted survivor function for a Cox proportional hazards
model with elastic net penalty from a cross-validated glmnet model.
}
\details{
This function makes it easier to use the results of cross-validation
to compute a survival curve.
}
\examples{
set.seed(2)
nobs <- 100; nvars <- 15
xvec <- rnorm(nobs * nvars)
x <- matrix(xvec, nrow = nobs)
beta <- rnorm(nvars / 3)
fx <- x[, seq(nvars / 3)] \%*\% beta / 3
ty <- rexp(nobs, exp(fx))
tcens <- rbinom(n = nobs, prob = 0.3, size = 1)
y <- survival::Surv(ty, tcens)
cvfit <- cv.glmnet(x, y, family = "cox")
# default: s = "lambda.1se"
survival::survfit(cvfit, x = x, y = y)

# s = "lambda.min"
survival::survfit(cvfit, s = "lambda.min", x = x, y = y)
}

\name{cv.glmnet}
\alias{cv.glmnet}
\title{Cross-validation for glmnet}
\description{Does k-fold cross-validation for glmnet, produces a plot,
  and returns a value for \code{lambda}}
\usage{
cv.glmnet(x, y, weights, offset, lambda, type.measure, nfolds, foldid, grouped, keep,
      parallel, ...)
}
\arguments{
  \item{x}{\code{x} matrix as in \code{glmnet}.}
  \item{y}{response \code{y} as in \code{glmnet}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{offset}{Offset vector (matrix) as in \code{glmnet}}
  \item{lambda}{Optional user-supplied lambda sequence; default is
    \code{NULL}, and \code{glmnet} chooses its own sequence}
  \item{nfolds}{number of folds - default is 10. Although \code{nfolds}
    can be as large as the sample size (leave-one-out CV), it is not
    recommended for large datasets. Smallest value allowable is \code{nfolds=3}}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing.}
  \item{type.measure}{loss to use for cross-validation. Currently five
    options, not all available for all models. The default is \code{type.measure="deviance"}, which uses
    squared-error for gaussian models (a.k.a \code{type.measure="mse"} there), deviance
    for logistic and poisson
    regression, and partial-likelihood for the Cox
    model. \code{type.measure="class"} applies to binomial and multinomial logistic regression only,
    and gives misclassification error.  \code{type.measure="auc"} is for
    two-class logistic regression only, and gives area under the ROC
    curve. \code{type.measure="mse"} or \code{type.measure="mae"} (mean absolute error)
    can be used by all models except the \code{"cox"}; they measure the
    deviation from the fitted mean to the response.} 
    \item{grouped}{This is an experimental argument, with default
      \code{TRUE}, and can be ignored by most users. For all models
      except the \code{"cox"}, this refers to computing \code{nfolds}
      separate statistics, and then using their mean and estimated
      standard error to describe the CV curve. If \code{grouped=FALSE},
      an error matrix is built up at the observation level from the predictions
      from the \code{nfold} fits, and then summarized (does not apply to
      \code{type.measure="auc"}). For the \code{"cox"} family,
      \code{grouped=TRUE} obtains the CV partial likelihood for the Kth
      fold by \emph{subtraction}; by subtracting the log partial
      likelihood evaluated on the full dataset from that evaluated on
      the  on the (K-1)/K dataset. This makes more efficient use of risk
      sets. With \code{grouped=FALSE} the log partial likelihood is
      computed only on the Kth fold}
    \item{keep}{If \code{keep=TRUE}, a \emph{prevalidated} array is
      returned containing fitted values for each observation and each
      value of \code{lambda}. This means these fits are computed with
      this observation and the rest of its fold omitted. The
      \code{folid} vector is also returned. Default is {keep=FALSE}}
    \item{parallel}{If \code{TRUE}, use parallel \code{foreach} to fit each fold.
      Must register parallel before hand, such as \code{doMC} or others.
    See the example below.}
 \item{\dots}{Other arguments that can be passed to \code{glmnet}}
}
\details{The function runs \code{glmnet} \code{nfolds}+1 times; the
  first to get the \code{lambda} sequence, and then the remainder to
  compute the fit with each of the folds omitted. The error is
  accumulated, and the average error and standard deviation over the
  folds is computed. 
  Note that \code{cv.glmnet} does NOT search for
  values for \code{alpha}. A specific value should be supplied, else
  \code{alpha=1} is assumed by default. If users would like to
  cross-validate \code{alpha} as well, they should call \code{cv.glmnet}
  with a pre-computed vector \code{foldid}, and then use this same fold vector
  in separate calls to \code{cv.glmnet} with different values of
  \code{alpha}. 
}
\value{an object of class \code{"cv.glmnet"} is returned, which is a
  list with the ingredients of the cross-validation fit.
  \item{lambda}{the values of \code{lambda} used in the fits.}
  \item{cvm}{The mean cross-validated error - a vector of length
    \code{length(lambda)}.}
  \item{cvsd}{estimate of standard error of \code{cvm}.}
  \item{cvup}{upper curve = \code{cvm+cvsd}.}
  \item{cvlo}{lower curve = \code{cvm-cvsd}.}
  \item{nzero}{number of non-zero coefficients at each \code{lambda}.}
  \item{name}{a text string indicating type of measure (for plotting
    purposes).}
 \item{glmnet.fit}{a fitted glmnet object for the full data.}
  \item{lambda.min}{value of \code{lambda} that gives minimum
    \code{cvm}.}
  \item{lambda.1se}{largest value of \code{lambda} such that error is
    within 1 standard error of the minimum.}
  \item{fit.preval}{if \code{keep=TRUE}, this is the array of
    prevalidated fits. Some entries can be \code{NA}, if that and
    subsequent values of \code{lambda} are not reached for that fold}
  \item{foldid}{if \code{keep=TRUE}, the fold assignments used}
}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
    Descent},   \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}\cr
  \emph{Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr
  Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011)
  \emph{Regularization Paths for Cox's Proportional Hazards Model via
    Coordinate Descent, Journal of Statistical Software, Vol. 39(5)
    1-13}\cr
  \url{http://www.jstatsoft.org/v39/i05/}
}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr 
  Noah Simon helped develop the 'coxnet' function.\cr
  Jeffrey Wong and B. Narasimhan helped with the parallel option\cr
Maintainer: Trevor Hastie \email{hastie@stanford.edu}}
\seealso{\code{glmnet} and \code{plot}, \code{predict}, and \code{coef} methods for \code{"cv.glmnet"} object.}
\examples{
set.seed(1010)
n=1000;p=100
nzc=trunc(p/10)
x=matrix(rnorm(n*p),n,p)
beta=rnorm(nzc)
fx= x[,seq(nzc)] \%*\% beta
eps=rnorm(n)*5
y=drop(fx+eps)
px=exp(fx)
px=px/(1+px)
ly=rbinom(n=length(px),prob=px,size=1)
set.seed(1011)
cvob1=cv.glmnet(x,y)
plot(cvob1)
coef(cvob1)
predict(cvob1,newx=x[1:5,], s="lambda.min")
title("Gaussian Family",line=2.5)
set.seed(1011)
cvob1a=cv.glmnet(x,y,type.measure="mae")
plot(cvob1a)
title("Gaussian Family",line=2.5)
set.seed(1011)
par(mfrow=c(2,2),mar=c(4.5,4.5,4,1))
cvob2=cv.glmnet(x,ly,family="binomial")
plot(cvob2)
title("Binomial Family",line=2.5)
frame()
set.seed(1011)
cvob3=cv.glmnet(x,ly,family="binomial",type.measure="class")
plot(cvob3)
title("Binomial Family",line=2.5)
set.seed(1011)
cvob3a=cv.glmnet(x,ly,family="binomial",type.measure="auc")
plot(cvob3a)
title("Binomial Family",line=2.5)
set.seed(1011)
mu=exp(fx/10)
y=rpois(n,mu)
cvob4=cv.glmnet(x,y,family="poisson")
plot(cvob4)
title("Poisson Family",line=2.5)
# Multinomial
n=500;p=30
nzc=trunc(p/10)
x=matrix(rnorm(n*p),n,p)
beta3=matrix(rnorm(30),10,3)
beta3=rbind(beta3,matrix(0,p-10,3))
f3=x\%*\% beta3
p3=exp(f3)
p3=p3/apply(p3,1,sum)
g3=rmult(p3)
set.seed(10101)
cvfit=cv.glmnet(x,g3,family="multinomial")
plot(cvfit)
title("Multinomial Family",line=2.5)
# Cox
beta=rnorm(nzc)
fx=x[,seq(nzc)]\%*\%beta/3
hx=exp(fx)
ty=rexp(n,hx)
tcens=rbinom(n=n,prob=.3,size=1)# censoring indicator
y=cbind(time=ty,status=1-tcens) # y=Surv(ty,1-tcens) with library(survival)
foldid=sample(rep(seq(10),length=n))
fit1_cv=cv.glmnet(x,y,family="cox",foldid=foldid)
plot(fit1_cv)
title("Cox Family",line=2.5)
\dontrun{
# Parallel
require(doMC)
registerDoMC(cores=4)
x = matrix(rnorm(1e5 * 100), 1e5, 100)
y = rnorm(1e5)
system.time(cv.glmnet(x,y))
system.time(cv.glmnet(x,y,parallel=TRUE))
}
}
\keyword{models}
\keyword{regression}


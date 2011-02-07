\name{glmnet}
\alias{glmnet}
\title{fit a GLM with lasso or elasticnet regularization}
\description{
    Fit a generalized linear model via penalized maximum likelihood.  The
  regularization path is computed for the lasso or elasticnet penalty at a grid
  of values for the regularization parameter lambda. Can deal with all
  shapes of data, including very large sparse data matrices. Fits
  linear, logistic and multinomial, poisson, and Cox regression models.
  }
\usage{
glmnet(x, y,
family=c("gaussian","binomial","poisson","multinomial","cox"), weights, offset=NULL,
alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs<nvars,0.01,0.0001), lambda=NULL,
standardize = TRUE, thresh = 1e-06,  dfmax = nvars + 1, pmax = min(dfmax * 1.2, nvars),
 exclude, penalty.factor = rep(1, nvars), maxit=100, HessianExact = FALSE, type.gaussian=ifelse(nvars<500,"covariance","naive"))
}

\arguments{
  \item{x}{input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse column format (class \code{"dgCMatrix"} as in package \code{Matrix}; not yet available for \code{family="cox"})}
  \item{y}{response variable. Quantitative for \code{family="gaussian"}
  or \code{family="poisson"} (non-negative counts). For
  \code{family="binomial"} should be either a factor with two levels, or
  a two-column matrix of counts or proportions. For
  \code{family="multinomial"}, can be a \code{nc>=2} level factor, or a
  matrix with \code{nc} columns of counts or proportions. For
  \code{family="cox"}, \code{y} should be a two-column matrix with
  columns named 'time' and 'status'. The latter is a binary variable,
  with '1' indicating death, and '0' indicating right censored. The
  function \code{Surv()} in package \pkg{survival} produces such a matrix}
  \item{family}{Response type (see above)}
  \item{weights}{observation weights. Can be total counts if responses are proportion matrices. Default is 1 for each observation}
  \item{offset}{A vector of length \code{nobs} that is included in the linear predictor (a \code{nobs x nc} matrix for the \code{"multinomial"} family). Useful for the \code{"poisson"} family (e.g. log of exposure time), or for refining a model by starting at a current fit. Default is \code{NULL}. If supplied, then values must also be supplied to the \code{predict} function.}
  \item{alpha}{The elasticnet mixing parameter, with
    \eqn{0<\alpha\le 1}. The penalty is defined
    as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1}
    is the lasso penalty.}
  \item{nlambda}{The number of \code{lambda} values - default is 100.}
  \item{lambda.min.ratio}{Smallest value for \code{lambda}, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero). The default depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs > nvars}, the default is \code{0.0001},
  close to zero.  If \code{nobs < nvars}, the default is \code{0.01}.
  A very small value of
  \code{lambda.min.ratio} will lead to a saturated fit in the \code{nobs <
  nvars} case. This is undefined for
  \code{"binomial"} and \code{"multinomial"} models, and \code{glmnet}
  will exit gracefully when the percentage deviance explained is almost
  1.} 
  \item{lambda}{A user supplied \code{lambda} sequence. Typical usage
    is to have the 
    program compute its own \code{lambda} sequence based on
    \code{nlambda} and \code{lambda.min.ratio}. Supplying a value of
    \code{lambda} overrides this. WARNING: use with care. Do not supply
  a single value for \code{lambda} (for predictions after CV use \code{predict()}
  instead).  Supply instead
    a decreasing sequence of \code{lambda} values. \code{glmnet} relies
  on its warms starts for speed, and its often faster to fit a whole
  path than compute a single fit.}
  \item{standardize}{Logical flag for variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. Default is is \code{standardize=TRUE}.
  If variables are in the same units already, you might not wish to standardize.}
  \item{thresh}{Convergence threshold for coordinate descent. Each inner
  coordinate-descent loop continues until the maximum change in the
  objective after any coefficient update is less than \code{thresh}
  times the null deviance. Defaults value is \code{1E-6}.}
  \item{dfmax}{Limit the maximum number of variables in the
    model. Useful for very large \code{nvars}, if a partial path is desired.}
  \item{pmax}{Limit the maximum number of variables ever to be nonzero}
  \item{exclude}{Indices of variables to be excluded from the
    model. Default is none. Equivalent to an infinite penalty factor
    (next item).}
  \item{penalty.factor}{Separate penalty factors can be applied to each
    coefficient. This is a number that multiplies \code{lambda} to allow
    differential shrinkage. Can be 0 for some variables, which implies
    no shrinkage, and that variable is always included in the
    model. Default is 1 for all variables (and implicitly infinity for
    variables listed in \code{exclude}).}
    \item{maxit}{Maximum number of outer-loop iterations for
      \code{"binomial"} or \code{"multinomial"} families. Default is 100.}
    \item{HessianExact}{Only applies to \code{"binomial"} or
	  \code{"multinomial"} families.
	  If \code{FALSE} (the default), an upper-bound approximation is made to the
    hessian, which is not recalculated at each outer loop.} 
  \item{type.gaussian}{Two algorithm types are supported for (only)
    \code{family="gaussian"}. The default when \code{nvar<500} is
  \code{type.gaussian="covariance"}, and saves all
    inner-products ever computed. This  can be much faster than
    \code{type.gaussian="naive"}, which loops through \code{nobs} every
  time an inner-product is computed. The latter can be far more efficient for \code{nvar >>
      nobs} situations, or when \code{nvar > 500}.}
}
\details{
  The sequence of models implied by \code{lambda} is fit by coordinate
  descent. For \code{family="gaussian"} this is the lasso sequence if
  \code{alpha=1}, else it is the elasticnet sequence.
 For the other families, this is a lasso or elasticnet regularization path
  for fitting the generalized linear regression
  paths, by maximizing the appropriate penalized log-likelihood (partial likelihood for the "cox" model). Sometimes the sequence is truncated before \code{nlambda}
  values of \code{lambda} have been used, because of instabilities in
  the inverse link functions near a saturated fit. \code{glmnet(...,family="binomial")}
  fits a traditional logistic regression model for the
  log-odds. \code{glmnet(...,family="multinomial")} fits a symmetric multinomial model, where
  each class is represented by a linear model (on the log-scale). The
  penalties take care of redundancies. A two-class \code{"multinomial"} model
  will produce the same fit as the corresponding \code{"binomial"} model,
  except the pair of coefficient matrices will be equal in magnitude and
  opposite in sign, and half the \code{"binomial"} values. 
  Note that the objective function for \code{"gaussian"} is \deqn{1/2 RSS/nobs + \lambda*penalty,} and for the other models it is  \deqn{-loglik/nobs + \lambda*penalty}
}
\value{
An object with S3 class \code{"glmnet","*" }, where \code{"*"} is
\code{"elnet"}, \code{"lognet"}, 
\code{"multnet"}, \code{"fishnet"} (poisson) or \code{"coxnet"}  for the five types of models.
  \item{call}{the call that produced this object}
  \item{a0}{Intercept sequence of length \code{length(lambda)}}
  \item{beta}{For \code{"elnet"} and \code{"lognet"} models, a \code{nvars x
      length(lambda)} matrix of coefficients, stored in sparse column
    format (\code{"dgCMatrix"}). For \code{"multnet"}, a list of \code{nc} such
    matrices, one for each class.}
  \item{lambda}{The actual sequence of \code{lambda} values used}
  \item{dev.ratio}{The fraction of (null) deviance explained (for \code{"elnet"}, this
      is the R-square). The deviance calculations incorporate weights if
  present in the model. The deviance is defined to be -2*(loglike_sat -
  loglike), where loglike_sat is the log-likelihood for the saturated
  model (a model with a free parameter per observation). Hence dev.fraction=1-dev/nulldev.}
    \item{nulldev}{Null deviance (per observation). This is defined to
  be  -2*(loglike_sat -loglike(Null)); The NULL model refers to the
  intercept model, except for the Cox, where it is the 0 model.}
  \item{df}{The number of nonzero coefficients for each value of
    \code{lambda}. For \code{"multnet"}, this is the number of variables
    with a nonzero coefficient for \emph{any} class.}
  \item{dfmat}{For \code{"multnet"} only. A matrix consisting of the
    number of nonzero coefficients per class}
  \item{dim}{dimension of coefficient matrix (ices)}
  \item{nobs}{number of observations}
  \item{npasses}{total passes over the data summed over all lambda
    values}
  \item{offset}{a logical variable indicating whether an offset was included in the model}
  \item{jerr}{error flag, for warnings and errors (largely for internal debugging).}
}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
    Descent, Journal of Statistical Software, Vol. 33, Issue 1, Feb 2010}\cr
  \url{http://www-stat.stanford.edu/~hastie/Papers/glmnet.pdf}\cr
  \url{http://www.jstatsoft.org/v33/i01/}
}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr 
Noah Simon helped develop the 'cox' feature.\cr
Maintainer: Trevor Hastie \email{hastie@stanford.edu}}
\seealso{\code{print}, \code{predict}, \code{coef} and \code{plot} methods, and the \code{cv.glmnet} function.}
\examples{
# Gaussian
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
fit1=glmnet(x,y)
print(fit1)
coef(fit1,s=0.01) # extract coefficients at a single value of lambda
predict(fit1,newx=x[1:10,],s=c(0.01,0.005)) # make predictions

#binomial
g2=sample(1:2,100,replace=TRUE)
fit2=glmnet(x,g2,family="binomial")

#multinomial
g4=sample(1:4,100,replace=TRUE)
fit3=glmnet(x,g4,family="multinomial")

#poisson
N=500; p=20
nzc=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
f = x[,seq(nzc)]\%*\%beta
mu=exp(f)
y=rpois(N,mu)
fit=glmnet(x,y,family="poisson")
plot(fit)
pfit = predict(fit,x,s=0.001,type="response")
plot(pfit,y)

#Cox
set.seed(10101)
N=1000;p=30
nzc=p/3
x=matrix(rnorm(N*p),N,p)
beta=rnorm(nzc)
fx=x[,seq(nzc)]\%*\%beta/3
hx=exp(fx)
ty=rexp(N,hx)
tcens=rbinom(n=N,prob=.3,size=1)# censoring indicator
y=cbind(time=ty,status=1-tcens) # y=Surv(ty,1-tcens) with library(survival)
fit=glmnet(x,y,family="cox")
plot(fit)

# Sparse
n=10000;p=200
nzc=trunc(p/10)
x=matrix(rnorm(n*p),n,p)
iz=sample(1:(n*p),size=n*p*.85,replace=FALSE)
x[iz]=0
sx=Matrix(x,sparse=TRUE)
beta=rnorm(nzc)
fx=x[,seq(nzc)]\%*\%beta
eps=rnorm(n)
y=fx+eps
px=exp(fx)
px=px/(1+px)
ly=rbinom(n=length(px),prob=px,size=1)
system.time(fit1<-glmnet(sx,y))
system.time(fit2n<-glmnet(x,y))
}
\keyword{models}
\keyword{regression}

 

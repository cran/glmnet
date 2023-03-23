#' Cross-validation for glmnet
#'
#' Does k-fold cross-validation for glmnet, produces a plot, and returns a
#' value for \code{lambda} (and \code{gamma} if \code{relax=TRUE})
#'
#' The function runs \code{glmnet} \code{nfolds}+1 times; the first to get the
#' \code{lambda} sequence, and then the remainder to compute the fit with each
#' of the folds omitted. The error is accumulated, and the average error and
#' standard deviation over the folds is computed.  Note that \code{cv.glmnet}
#' does NOT search for values for \code{alpha}. A specific value should be
#' supplied, else \code{alpha=1} is assumed by default. If users would like to
#' cross-validate \code{alpha} as well, they should call \code{cv.glmnet} with
#' a pre-computed vector \code{foldid}, and then use this same fold vector in
#' separate calls to \code{cv.glmnet} with different values of \code{alpha}.
#' Note also that the results of \code{cv.glmnet} are random, since the folds
#' are selected at random. Users can reduce this randomness by running
#' \code{cv.glmnet} many times, and averaging the error curves.
#'
#' If \code{relax=TRUE} then the values of \code{gamma} are used to mix the
#' fits. If \eqn{\eta} is the fit for lasso/elastic net, and \eqn{\eta_R} is
#' the relaxed fit (with unpenalized coefficients), then a relaxed fit mixed by
#' \eqn{\gamma} is \deqn{\eta(\gamma)=(1-\gamma)\eta_R+\gamma\eta.} There is
#' practically no extra cost for having a lot of values for \code{gamma}.
#' However, 5 seems sufficient for most purposes. CV then selects both
#' \code{gamma} and \code{lambda}.
#'
#' @param x \code{x} matrix as in \code{glmnet}.
#' @param y response \code{y} as in \code{glmnet}.
#' @param weights Observation weights; defaults to 1 per observation
#' @param offset Offset vector (matrix) as in \code{glmnet}
#' @param lambda Optional user-supplied lambda sequence; default is
#' \code{NULL}, and \code{glmnet} chooses its own sequence. Note that this is done
#' for the full model (master sequence), and separately for each fold.
#' The fits are then alligned using the master sequence (see the \code{allignment}
#' argument for additional details). Adapting \code{lambda} for each fold
#' leads to better convergence. When \code{lambda} is supplied, the same sequence
#' is used everywhere, but in some GLMs can lead to convergence issues.
#' @param type.measure loss to use for cross-validation. Currently five
#' options, not all available for all models. The default is
#' \code{type.measure="deviance"}, which uses squared-error for gaussian models
#' (a.k.a \code{type.measure="mse"} there), deviance for logistic and poisson
#' regression, and partial-likelihood for the Cox model.
#' \code{type.measure="class"} applies to binomial and multinomial logistic
#' regression only, and gives misclassification error.
#' \code{type.measure="auc"} is for two-class logistic regression only, and
#' gives area under the ROC curve. \code{type.measure="mse"} or
#' \code{type.measure="mae"} (mean absolute error) can be used by all models
#' except the \code{"cox"}; they measure the deviation from the fitted mean to
#' the response.
#' \code{type.measure="C"} is Harrel's concordance measure, only available for \code{cox} models.
#' @param nfolds number of folds - default is 10. Although \code{nfolds} can be
#' as large as the sample size (leave-one-out CV), it is not recommended for
#' large datasets. Smallest value allowable is \code{nfolds=3}
#' @param foldid an optional vector of values between 1 and \code{nfolds}
#' identifying what fold each observation is in. If supplied, \code{nfolds} can
#' be missing.
#' @param alignment This is an experimental argument, designed to fix the
#' problems users were having with CV, with possible values \code{"lambda"}
#' (the default) else \code{"fraction"}. With \code{"lambda"} the \code{lambda}
#' values from the master fit (on all the data) are used to line up the
#' predictions from each of the folds. In some cases this can give strange
#' values, since the effective \code{lambda} values in each fold could be quite
#' different. With \code{"fraction"} we line up the predictions in each fold
#' according to the fraction of progress along the regularization. If in the
#' call a \code{lambda} argument is also provided, \code{alignment="fraction"}
#' is ignored (with a warning).
#' @param grouped This is an experimental argument, with default \code{TRUE},
#' and can be ignored by most users. For all models except the \code{"cox"},
#' this refers to computing \code{nfolds} separate statistics, and then using
#' their mean and estimated standard error to describe the CV curve. If
#' \code{grouped=FALSE}, an error matrix is built up at the observation level
#' from the predictions from the \code{nfolds} fits, and then summarized (does
#' not apply to \code{type.measure="auc"}). For the \code{"cox"} family,
#' \code{grouped=TRUE} obtains the CV partial likelihood for the Kth fold by
#' \emph{subtraction}; by subtracting the log partial likelihood evaluated on
#' the full dataset from that evaluated on the on the (K-1)/K dataset. This
#' makes more efficient use of risk sets. With \code{grouped=FALSE} the log
#' partial likelihood is computed only on the Kth fold
#' @param keep If \code{keep=TRUE}, a \emph{prevalidated} array is returned
#' containing fitted values for each observation and each value of
#' \code{lambda}. This means these fits are computed with this observation and
#' the rest of its fold omitted. The \code{foldid} vector is also returned.
#' Default is keep=FALSE.  If \code{relax=TRUE}, then a list of such arrays is
#' returned, one for each value of 'gamma'. Note: if the value 'gamma=1' is
#' omitted, this case is included in the list since it corresponds to the
#' original 'glmnet' fit.
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to fit each
#' fold.  Must register parallel before hand, such as \code{doMC} or others.
#' See the example below.
#' @param gamma The values of the parameter for mixing the relaxed fit with the
#' regularized fit, between 0 and 1; default is \code{gamma = c(0, 0.25, 0.5,
#' 0.75, 1)}
#' @param relax If \code{TRUE}, then CV is done with respect to the mixing
#' parameter \code{gamma} as well as \code{lambda}. Default is
#' \code{relax=FALSE}
#' @param trace.it If \code{trace.it=1}, then progress bars are displayed;
#' useful for big models that take a long time to fit. Limited tracing if
#' \code{parallel=TRUE}
#' @param \dots Other arguments that can be passed to \code{glmnet}
#' @return an object of class \code{"cv.glmnet"} is returned, which is a list
#' with the ingredients of the cross-validation fit.  If the object was created
#' with \code{relax=TRUE} then this class has a prefix class of
#' \code{"cv.relaxed"}.  \item{lambda}{the values of \code{lambda} used in the
#' fits.} \item{cvm}{The mean cross-validated error - a vector of length
#' \code{length(lambda)}.} \item{cvsd}{estimate of standard error of
#' \code{cvm}.} \item{cvup}{upper curve = \code{cvm+cvsd}.} \item{cvlo}{lower
#' curve = \code{cvm-cvsd}.} \item{nzero}{number of non-zero coefficients at
#' each \code{lambda}.} \item{name}{a text string indicating type of measure
#' (for plotting purposes).} \item{glmnet.fit}{a fitted glmnet object for the
#' full data.} \item{lambda.min}{value of \code{lambda} that gives minimum
#' \code{cvm}.} \item{lambda.1se}{largest value of \code{lambda} such that
#' error is within 1 standard error of the minimum.} \item{fit.preval}{if
#' \code{keep=TRUE}, this is the array of prevalidated fits. Some entries can
#' be \code{NA}, if that and subsequent values of \code{lambda} are not reached
#' for that fold} \item{foldid}{if \code{keep=TRUE}, the fold assignments used}
#' \item{index}{a one column matrix with the indices of \code{lambda.min} and \code{lambda.1se} in the sequence of coefficients, fits etc.}
#' \item{relaxed}{if \code{relax=TRUE}, this additional item has the CV info
#' for each of the mixed fits. In particular it also selects \code{lambda,
#' gamma} pairs corresponding to the 1se rule, as well as the minimum error. It also has a component \code{index}, a two-column matrix  which contains the \code{lambda} and \code{gamma} indices corresponding to the "min" and "1se" solutions.}
#' @author Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr Noah Simon
#' helped develop the 'coxnet' function.\cr Jeffrey Wong and B. Narasimhan
#' helped with the parallel option\cr Maintainer: Trevor Hastie
#' \email{hastie@@stanford.edu}
#' @seealso \code{glmnet} and \code{plot}, \code{predict}, and \code{coef}
#' methods for \code{"cv.glmnet"} and \code{"cv.relaxed"} objects.
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#' \emph{Regularization Paths for Generalized Linear Models via Coordinate
#' Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22},
#' \doi{10.18637/jss.v033.i01}.\cr
#' Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011)
#' \emph{Regularization Paths for Cox's Proportional
#' Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol.
#' 39(5), 1-13},
#' \doi{10.18637/jss.v039.i05}.
#' @keywords models regression
#' @examples
#'
#' set.seed(1010)
#' n = 1000
#' p = 100
#' nzc = trunc(p/10)
#' x = matrix(rnorm(n * p), n, p)
#' beta = rnorm(nzc)
#' fx = x[, seq(nzc)] %*% beta
#' eps = rnorm(n) * 5
#' y = drop(fx + eps)
#' px = exp(fx)
#' px = px/(1 + px)
#' ly = rbinom(n = length(px), prob = px, size = 1)
#' set.seed(1011)
#' cvob1 = cv.glmnet(x, y)
#' plot(cvob1)
#' coef(cvob1)
#' predict(cvob1, newx = x[1:5, ], s = "lambda.min")
#' title("Gaussian Family", line = 2.5)
#' set.seed(1011)
#' cvob1a = cv.glmnet(x, y, type.measure = "mae")
#' plot(cvob1a)
#' title("Gaussian Family", line = 2.5)
#' set.seed(1011)
#' par(mfrow = c(2, 2), mar = c(4.5, 4.5, 4, 1))
#' cvob2 = cv.glmnet(x, ly, family = "binomial")
#' plot(cvob2)
#' title("Binomial Family", line = 2.5)
#' frame()
#' set.seed(1011)
#' cvob3 = cv.glmnet(x, ly, family = "binomial", type.measure = "class")
#' plot(cvob3)
#' title("Binomial Family", line = 2.5)
#' \dontrun{
#' cvob1r = cv.glmnet(x, y, relax = TRUE)
#' plot(cvob1r)
#' predict(cvob1r, newx = x[, 1:5])
#' set.seed(1011)
#' cvob3a = cv.glmnet(x, ly, family = "binomial", type.measure = "auc")
#' plot(cvob3a)
#' title("Binomial Family", line = 2.5)
#' set.seed(1011)
#' mu = exp(fx/10)
#' y = rpois(n, mu)
#' cvob4 = cv.glmnet(x, y, family = "poisson")
#' plot(cvob4)
#' title("Poisson Family", line = 2.5)
#'
#' # Multinomial
#' n = 500
#' p = 30
#' nzc = trunc(p/10)
#' x = matrix(rnorm(n * p), n, p)
#' beta3 = matrix(rnorm(30), 10, 3)
#' beta3 = rbind(beta3, matrix(0, p - 10, 3))
#' f3 = x %*% beta3
#' p3 = exp(f3)
#' p3 = p3/apply(p3, 1, sum)
#' g3 = glmnet:::rmult(p3)
#' set.seed(10101)
#' cvfit = cv.glmnet(x, g3, family = "multinomial")
#' plot(cvfit)
#' title("Multinomial Family", line = 2.5)
#' # Cox
#' beta = rnorm(nzc)
#' fx = x[, seq(nzc)] %*% beta/3
#' hx = exp(fx)
#' ty = rexp(n, hx)
#' tcens = rbinom(n = n, prob = 0.3, size = 1)  # censoring indicator
#' y = cbind(time = ty, status = 1 - tcens)  # y=Surv(ty,1-tcens) with library(survival)
#' foldid = sample(rep(seq(10), length = n))
#' fit1_cv = cv.glmnet(x, y, family = "cox", foldid = foldid)
#' plot(fit1_cv)
#' title("Cox Family", line = 2.5)
#' # Parallel
#' require(doMC)
#' registerDoMC(cores = 4)
#' x = matrix(rnorm(1e+05 * 100), 1e+05, 100)
#' y = rnorm(1e+05)
#' system.time(cv.glmnet(x, y))
#' system.time(cv.glmnet(x, y, parallel = TRUE))
#' }
#'
#' @export cv.glmnet
cv.glmnet <-
  function (x, y, weights=NULL, offset = NULL, lambda = NULL, type.measure = c("default","mse",
                                                           "deviance", "class", "auc", "mae","C"), nfolds = 10, foldid=NULL,  alignment=c("lambda","fraction"),grouped = TRUE, keep = FALSE, parallel = FALSE, gamma=c(0,.25,.5,.75,1),relax=FALSE,trace.it=0, ...)
{
    type.measure = match.arg(type.measure)
    alignment=match.arg(alignment)
  if (!is.null(lambda) && length(lambda) < 2)
      stop("Need more than one value of lambda for cv.glmnet")
  if (!is.null(lambda) && alignment=="fraction"){
      warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
      alignment="lambda"
  }
   N = nrow(x)
 if (is.null(weights))
    weights = rep(1, N)
  else weights = as.double(weights)
  y = drop(y)
  cv.call=glmnet.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", "grouped",
    "keep"), names(glmnet.call), FALSE)
  if (any(which))
    glmnet.call = glmnet.call[-which]
    glmnet.call[[1]] = as.name("glmnet")
    if(glmnet.control()$itrace)trace.it=1
    else{
        if(trace.it){
            glmnet.control(itrace=1)
            on.exit(glmnet.control(itrace=0))
            }
        }
   if (is.null(foldid))
    foldid = sample(rep(seq(nfolds), length = N))
  else nfolds = max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
### Now we switch depending on relax
    if(relax)
        cv.relaxed.raw(x,y,weights,offset,lambda,type.measure,nfolds,foldid,
                       alignment,grouped,keep,parallel,trace.it,glmnet.call,cv.call,gamma,...)
    else
        cv.glmnet.raw(x,y,weights,offset,lambda,type.measure,nfolds,foldid,
                       alignment,grouped,keep,parallel,trace.it,glmnet.call,cv.call,...)
 }


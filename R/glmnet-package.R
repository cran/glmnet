

#' Simulated data for the glmnet vignette
#'
#' Simple simulated data, used to demonstrate the features of glmnet
#'
#' These datasets are artificial, and are used to test out some of the
#' features of glmnet.
#' @name beta_CVX
#' @aliases x y beta_CVX
#' @format Data objects used to demonstrate features in the glmnet vignette
#' @keywords datasets
#' @useDynLib glmnet
#' @import methods
#' @import Matrix
#' @import foreach
#' @importFrom utils packageDescription
#' @importFrom graphics abline axis matplot points segments text par plot
#' @importFrom stats  approx  coef  median  predict  rmultinom  runif  weighted.mean family rnorm gaussian binomial glm glm.fit
#' @importFrom survival concordance Surv is.Surv
#' @importFrom grDevices rainbow
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#'
#' data(QuickStartExample)
#' x <- QuickStartExample$x; y <- QuickStartExample$y
#' glmnet(x, y)
#'
NULL


#' Internal glmnet functions
#'
#' @description
#' These are not intended for use by users. \code{lambda.interp} does linear
#' interpolation of the lambdas to obtain a prediction at a new point s.
#' \code{glmnet_softmax} does the classification for multinomial models.
#' \code{nonzeroCoef} determines in an efficient manner which variables are
#' nonzero in each fit. \code{jerr} prints out error messages from the C++ routines.
#' \code{plotCoef} is called by the \code{plot} method for \code{glmnet}
#' objects. \code{check_dots} is used in \code{coef} and \code{predict} with
#' argument \code{exact=TRUE}, to make sure user supplies original data used to
#' fit the \code{"glmnet"} object.
#'
#' @name glmnet-internal
#' @aliases auc assess.coxnet auc.mat cvtype cvstats
#' cvcompute getcoef getcoef.multinomial fix.lam error.bars
#' getmin elnet mrelnet lognet fishnet coefnorm coxnet cv.lognet cv.elnet
#' cv.multnet cv.mrelnet cv.coxnet cv.fishnet cv.glmnet.raw cv.relaxed.raw
#' blend.relaxed checkgamma.relax buildPredmat buildPredmat.mrelnetlist
#' buildPredmat.multnetlist buildPredmat.lognetlist buildPredmat.array
#' buildPredmat.coxnetlist buildPredmat.default lambda.interp nonzeroCoef
#' glmnet_softmax getOptcv.glmnet getOptcv.relaxed jerr jerr.elnet jerr.lognet
#' jerr.fishnet jerr.coxnet jerr.mrelnet plotCoef zeromat na.mean check_dots
#' na_sparse_fix prepareX
#'
#' @author Trevor Hastie
#' @keywords internal
NULL





#' Elastic net model paths for some generalized linear models
#'
#' This package fits lasso and elastic-net model paths for regression, logistic
#' and multinomial regression using coordinate descent. The algorithm is
#' extremely fast, and exploits sparsity in the input x matrix where it exists.
#' A variety of predictions can be made from the fitted models.
#'
#' \tabular{ll}{ Package: \tab glmnet\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2008-05-14\cr License: \tab What license is it under?\cr }
#' Very simple to use. Accepts \code{x,y} data for regression models, and
#' produces the regularization path over a grid of values for the tuning
#' parameter \code{lambda}. Only 5 functions: \code{glmnet}\cr
#' \code{predict.glmnet}\cr \code{plot.glmnet}\cr \code{print.glmnet}\cr
#' \code{coef.glmnet}
#'
#' @name glmnet-package
#' @docType package
#' @author Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr Maintainer:
#' Trevor Hastie <hastie@@stanford.edu>
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#' \emph{Regularization Paths for Generalized Linear Models via Coordinate
#' Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22},
#' \doi{10.18637/jss.v033.i01}.\cr
#' Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011)
#' \emph{Regularization Paths for Cox's Proportional
#' Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol.
#' 39(5), 1-13},
#' \doi{10.18637/jss.v039.i05}.\cr
#' Tibshirani,Robert, Bien, J., Friedman, J., Hastie, T.,Simon, N.,Taylor, J. and
#' Tibshirani, Ryan. (2012) \emph{Strong Rules for Discarding Predictors in
#' Lasso-type Problems, JRSSB, Vol. 74(2), 245-266},
#' \url{https://arxiv.org/abs/1011.2234}.\cr
#' Hastie, T., Tibshirani, Robert and Tibshirani, Ryan (2020) \emph{Best Subset,
#' Forward Stepwise or Lasso? Analysis and Recommendations Based on Extensive Comparisons,
#' Statist. Sc. Vol. 35(4), 579-592},
#' \url{https://arxiv.org/abs/1707.08692}.\cr
#' Glmnet webpage with four vignettes: \url{https://glmnet.stanford.edu}.
#' @keywords models regression package
#' @examples
#'
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' g2 = sample(1:2, 100, replace = TRUE)
#' g4 = sample(1:4, 100, replace = TRUE)
#' fit1 = glmnet(x, y)
#' predict(fit1, newx = x[1:5, ], s = c(0.01, 0.005))
#' predict(fit1, type = "coef")
#' plot(fit1, xvar = "lambda")
#' fit2 = glmnet(x, g2, family = "binomial")
#' predict(fit2, type = "response", newx = x[2:5, ])
#' predict(fit2, type = "nonzero")
#' fit3 = glmnet(x, g4, family = "multinomial")
#' predict(fit3, newx = x[1:3, ], type = "response", s = 0.01)
#'
NULL




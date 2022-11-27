#' Check if glmnet should call cox.path
#'
#' Helper function to check if glmnet() should call cox.path().
#'
#' For \code{family="cox"}, we only call the original coxnet() function if
#' (i) x is not sparse, (ii) y is right-censored data, and (iii) we are
#' not fitting a stratified Cox model. This function also throws an error
#' if y has a "strata" attribute but is not of type "stratifySurv".
#'
#' @param x Design matrix.
#' @param y Response variable.
#'
#' @return TRUE if cox.path() should be called, FALSE otherwise.
use.cox.path <- function(x, y) {
  y <- response.coxnet(y)
  use_cox_path <- TRUE
  # We only return FALSE if:
  # x is not sparse AND y is right-censored data AND no strata
  # (strata variable being all equal counts as no strata)
  if (!inherits(x, "sparseMatrix") && ncol(y) == 2) {
    if (!("strata" %in% names(attributes(y))) ||
        length(unique(attr(y, "strata"))) == 1)
      use_cox_path <- FALSE
  }
  # if strata provided in y but y not of class stratifySurv, throw error
  if ("strata" %in% names(attributes(y)) && !inherits(y, "stratifySurv"))
    stop(paste0("For fitting stratified Cox models, y must be of class ",
                "stratifySurv, see ?stratifySurv for more details"))

  return(use_cox_path)
}

#' Fit a Cox regression model with elastic net regularization for a path of
#' lambda values
#'
#' Fit a Cox regression model via penalized maximum likelihood for a path of
#' lambda values. Can deal with (start, stop] data and strata, as well as
#' sparse design matrices.
#'
#' Sometimes the sequence is truncated before \code{nlambda} values of lambda
#' have been used. This happens when \code{cox.path} detects that the
#' decrease in deviance is marginal (i.e. we are near a saturated fit).
#'
#' @param x See glmnet help file
#' @param y Survival response variable, must be a \code{Surv} or
#' \code{stratifySurv} object.
#' @param weights See glmnet help file
#' @param offset See glmnet help file
#' @param alpha See glmnet help file
#' @param nlambda See glmnet help file
#' @param lambda.min.ratio See glmnet help file
#' @param lambda See glmnet help file
#' @param standardize See glmnet help file
#' @param thresh Convergence threshold for coordinate descent. Each inner
#' coordinate-descent loop continues until the maximum change in the objective
#' after any coefficient update is less than thresh times the null deviance.
#' Default value is \code{1e-10}.
#' @param exclude See glmnet help file
#' @param penalty.factor See glmnet help file
#' @param lower.limits See glmnet help file
#' @param upper.limits See glmnet help file
#' @param maxit See glmnet help file
#' @param trace.it Controls how much information is printed to screen. Default is
#' \code{trace.it=0} (no information printed). If \code{trace.it=1}, a progress
#' bar is displayed. If \code{trace.it=2}, some information about the fitting
#' procedure is printed to the console as the model is being fitted.
#' @param ... Other arguments passed from glmnet (not used right now).
#'
#' @return An object of class "coxnet" and "glmnet".
#' \item{a0}{Intercept value, \code{NULL} for "cox" family.}
#' \item{beta}{A \code{nvars x length(lambda)} matrix of coefficients, stored in
#' sparse matrix format.}
#' \item{df}{The number of nonzero coefficients for each value of lambda.}
#' \item{dim}{Dimension of coefficient matrix.}
#' \item{lambda}{The actual sequence of lambda values used. When alpha=0, the
#' largest lambda reported does not quite give the zero coefficients reported
#' (lambda=inf would in principle). Instead, the largest lambda for alpha=0.001
#' is used, and the sequence of lambda values is derived from this.}
#' \item{dev.ratio}{The fraction of (null) deviance explained. The deviance
#' calculations incorporate weights if present in the model. The deviance is
#' defined to be 2*(loglike_sat - loglike), where loglike_sat is the log-likelihood
#' for the saturated model (a model with a free parameter per observation).
#' Hence dev.ratio=1-dev/nulldev.}
#' \item{nulldev}{Null deviance (per observation). This is defined to be
#' 2*(loglike_sat -loglike(Null)). The null model refers to the 0 model.}
#' \item{npasses}{Total passes over the data summed over all lambda values.}
#' \item{jerr}{Error flag, for warnings and errors (largely for internal
#' debugging).}
#' \item{offset}{A logical variable indicating whether an offset was included
#' in the model.}
#' \item{call}{The call that produced this object.}
#' \item{nobs}{Number of observations.}
#'
#' @examples
#' set.seed(2)
#' nobs <- 100; nvars <- 15
#' xvec <- rnorm(nobs * nvars)
#' xvec[sample.int(nobs * nvars, size = 0.4 * nobs * nvars)] <- 0
#' x <- matrix(xvec, nrow = nobs)
#' beta <- rnorm(nvars / 3)
#' fx <- x[, seq(nvars / 3)] %*% beta / 3
#' ty <- rexp(nobs, exp(fx))
#' tcens <- rbinom(n = nobs, prob = 0.3, size = 1)
#' jsurv <- survival::Surv(ty, tcens)
#' fit1 <- glmnet:::cox.path(x, jsurv)
#'
#' # works with sparse x matrix
#' x_sparse <- Matrix::Matrix(x, sparse = TRUE)
#' fit2 <- glmnet:::cox.path(x_sparse, jsurv)
#'
#' # example with (start, stop] data
#' set.seed(2)
#' start_time <- runif(100, min = 0, max = 5)
#' stop_time <- start_time + runif(100, min = 0.1, max = 3)
#' status <- rbinom(n = nobs, prob = 0.3, size = 1)
#' jsurv_ss <- survival::Surv(start_time, stop_time, status)
#' fit3 <- glmnet:::cox.path(x, jsurv_ss)
#'
#' # example with strata
#' jsurv_ss2 <- stratifySurv(jsurv_ss, rep(1:2, each = 50))
#' fit4 <- glmnet:::cox.path(x, jsurv_ss2)
cox.path <- function(x, y, weights=NULL, offset=NULL,
                     alpha=1.0, nlambda=100,
                     lambda.min.ratio=ifelse(nobs<nvars, 1e-2, 1e-4),
                     lambda=NULL, standardize=TRUE,
                     thresh=1e-10, exclude=NULL, penalty.factor=rep(1,nvars),
                     lower.limits=-Inf, upper.limits=Inf, maxit=100000,
                     trace.it=0, ...) {
  ### Prepare all the generic arguments (mimicking top-level glmnet() call)
  if (alpha > 1) {
    warning("alpha > 1; set to 1")
    alpha = 1
  } else if (alpha < 0) {
    warning("alpha < 0; set to 0")
    alpha = 0
  }
  alpha = as.double(alpha)

  this.call <- match.call()

  np = dim(x)
  if (is.null(np) || (np[2] <= 1)) stop("x should be a matrix with 2 or more columns")
  nobs = as.integer(np[1]); nvars = as.integer(np[2])

  # get feature variable names
  vnames <- colnames(x)
  if(is.null(vnames)) vnames <- paste("V",seq(nvars),sep="")

  # check weights
  if(is.null(weights)) weights = rep(1,nobs)
  else if (length(weights) != nobs)
    stop(paste("Number of elements in weights (",length(weights),
               ") not equal to the number of rows of x (",nobs,")",sep=""))
  weights <- as.double(weights)

  # check that response y is a Surv object of the correct length
  y <- response.coxnet(y)
  if (nrow(y) != nobs) stop(paste0("number of observations in y (" , nrow(y),
                                   ") not equal to the number of rows of x (",
                                   nobs, ")"))

  # check offset option
  is.offset <- !(is.null(offset))
  if (is.offset == FALSE) {
    offset <- rep(0, times = nrow(y))
  }

  # check and standardize penalty factors (to sum to nvars)
  if(any(penalty.factor == Inf)) {
    exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
    exclude = sort(unique(exclude))
  }

    ## Compute weighted mean and variance of columns of x, sensitive to sparse matrix
    ## needed to detect constant columns below, and later if standarization
    meansd <- weighted_mean_sd(x, weights)

    ## look for constant variables, and if any, then add to exclude
    const_vars <- meansd$sd == 0
    nzvar <- setdiff(which(!const_vars), exclude)
    # if all the non-excluded variables have zero variance, throw error
    if (length(nzvar) == 0) stop("All used predictors have zero variance")

    ## if any constant vars, add to exclude
    if(any(const_vars)) {
        exclude <- sort(unique(c(which(const_vars),exclude)))
        meansd$sd[const_vars] <- 1.0 ## we divide later, and do not want bad numbers
    }
    if(length(exclude) > 0) {
        jd = match(exclude, seq(nvars), 0)
        if(!all(jd > 0)) stop ("Some excluded variables out of range")
        penalty.factor[jd] = 1 # ow can change lambda sequence
    }
    # check and standardize penalty factors (to sum to nvars)
    vp = pmax(0, penalty.factor)
    if (max(vp) <= 0) stop("All penalty factors are <= 0")
    vp = as.double(vp * nvars / sum(vp))

  ### check on limits
  control <- glmnet.control()
  if (thresh >= control$epsnr)
    warning("thresh should be smaller than glmnet.control()$epsnr",
            call. = FALSE)

  if(any(lower.limits > 0)){ stop("Lower limits should be non-positive") }
  if(any(upper.limits < 0)){ stop("Upper limits should be non-negative") }
  lower.limits[lower.limits == -Inf] = -control$big
  upper.limits[upper.limits == Inf] = control$big
  if (length(lower.limits) < nvars) {
    if(length(lower.limits) == 1) lower.limits = rep(lower.limits, nvars) else
      stop("Require length 1 or nvars lower.limits")
  } else lower.limits = lower.limits[seq(nvars)]
  if (length(upper.limits) < nvars) {
    if(length(upper.limits) == 1) upper.limits = rep(upper.limits, nvars) else
      stop("Require length 1 or nvars upper.limits")
  } else upper.limits = upper.limits[seq(nvars)]

  if (any(lower.limits == 0) || any(upper.limits == 0)) {
    ###Bounds of zero can mess with the lambda sequence and fdev;
    ###ie nothing happens and if fdev is not zero, the path can stop
    fdev <- glmnet.control()$fdev
    if(fdev!= 0) {
      glmnet.control(fdev = 0)
      on.exit(glmnet.control(fdev = fdev))
    }
  }
  ### end check on limits
  ### end preparation of generic arguments

  # standardize x if necessary
    xm <- rep(0.0, times = nvars)
    if (standardize) {
        xs <- meansd$sd
    } else {
        xs <- rep(1.0, times = nvars)
    }
    if (!inherits(x, "sparseMatrix")) {
    x <- scale(x,FALSE,xs)
  } else {
    attr(x, "xm") <- xm
    attr(x, "xs") <- xs
  }
  lower.limits <- lower.limits * xs
  upper.limits <- upper.limits * xs

  if (!("strata" %in% names(attributes(y))))
    y <- stratifySurv(y, rep(1, nobs))

  # Pre-compute and cache some important information: ordering by stop time
  # (ascending, deaths before censored), and for (start, stop] data: ordering
  # by start time and some match information.
  # Information is computed at the strata level
  if (ncol(y) == 2) {
    stop_o <- numeric(nobs)
    for (i in unique(attr(y, "strata"))) {
      ii <- which(attr(y, "strata") == i)
      stop_o[ii] <- order(y[ii, "time"], y[ii, "status"],
                          decreasing = c(FALSE, TRUE))
    }
    attr(y, "stop_time") <- stop_o
  } else {
    stop_o   <- numeric(nobs)
    start_o  <- numeric(nobs)
    ss_match <- numeric(nobs)
    for (i in unique(attr(y, "strata"))) {
      ii <- which(attr(y, "strata") == i)
      stop_o[ii]  <- order(y[ii, "stop"], y[ii, "status"],
                           decreasing = c(FALSE, TRUE))
      start_o[ii] <- order(y[ii, "start"], decreasing = c(FALSE))
      ss_match[ii] <- match(start_o[ii], stop_o[ii])
    }
    attr(y, "stop_time")  <- stop_o
    attr(y, "start_time") <- start_o
    attr(y, "ss_match")   <- ss_match
  }

  # compute null deviance
  # currently using std.weights = FALSE in order to match glmnet output
  nulldev <- coxnet.deviance(y = y, offset = offset, weights = weights,
                             std.weights = FALSE)

  # compute lambda_max and lambda values
  nlam = as.integer(nlambda)
  user_lambda = FALSE   # did user provide their own lambda values?
  if (is.null(lambda)) {
    if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
    lambda_max <- get_cox_lambda_max(x, y, alpha, weights, offset, exclude, vp)
    ulam <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                    length.out = nlam))
  } else {  # user provided lambda values
    user_lambda = TRUE
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }

  # start progress bar
  if (trace.it == 1) pb  <- createPB(min = 0, max = nlam, style = 3)

  beta <- matrix(0, nrow = nvars, ncol = nlam)
  dev.ratio <- rep(NA, length = nlam)
  fit <- NULL
  mnl <- min(nlam, control$mnlam)
  for (k in 1:nlam) {
    # get the correct lambda value to fit
    if (k > 1) {
      cur_lambda <- ulam[k]
    } else {
      cur_lambda <- ifelse(user_lambda, ulam[k], control$big)
    }

    if (trace.it == 2) cat("Fitting lambda index", k, ":", ulam[k], fill = TRUE)
    fit <- cox.fit(x, y, weights / sum(weights),
                   lambda = cur_lambda, alpha = alpha,
                   offset = offset, thresh = thresh, maxit = maxit,
                   penalty.factor = vp, exclude = exclude,
                   lower.limits = lower.limits, upper.limits = upper.limits,
                   warm = fit, from.cox.path = TRUE, save.fit = TRUE,
                   trace.it = trace.it)
    if (trace.it == 1) utils::setTxtProgressBar(pb, k)
    # if error code non-zero, a non-fatal error must have occurred
    # print warning, ignore this lambda value and return result
    # for all previous lambda values
    if (fit$jerr != 0) {
      errmsg <- jerr.glmnetfit(fit$jerr, maxit, k)
      warning(errmsg$msg, call. = FALSE)
      k <- k - 1
      break
    }
    beta[, k] <- as.vector(fit$beta)
    dev.ratio[k] <- fit$dev.ratio

    # early stopping if dev.ratio almost 1 or no improvement
    if (k >= mnl && user_lambda == FALSE) {
      if (dev.ratio[k] > control$devmax * 0.99 / 0.999) break
      if (k > 1 && dev.ratio[k] - dev.ratio[k - mnl + 1] <
          control$fdev * 100 * dev.ratio[k]) break
    }
  }
  if (trace.it == 1) {
    utils::setTxtProgressBar(pb, nlam)
    cat("", fill = TRUE)
  }

  # truncate beta, dev.ratio, lambda if necessary
  if (k < nlam) {
    beta <- beta[, 1:k, drop = FALSE]
    dev.ratio <- dev.ratio[1:k]
    ulam <- ulam[1:k]
  }

  # return coefficients to original scale (because of x standardization)
  beta <- beta / xs

  # output
  stepnames <- paste0("s", 0:(length(ulam) - 1))
  out <- list(a0 = NULL)
  out$beta <- Matrix::Matrix(beta, sparse = TRUE,
                             dimnames = list(vnames, stepnames))
  out$df <- as.vector(colSums(abs(beta) > 0))  # as.vector to remove names
  out$dim <- dim(beta)
  out$lambda <- ulam
  out$dev.ratio <- dev.ratio
  out$nulldev <- nulldev
  out$npasses <- fit$npasses
  out$jerr <- fit$jerr
  out$offset <- is.offset
  out$call <- this.call
  out$nobs <- nobs
  class(out) <- c("coxnet", "glmnet")
  out
}

#' Fit a Cox regression model with elastic net regularization for a single
#' value of lambda
#'
#' Fit a Cox regression model via penalized maximum likelihood for a single
#' value of lambda. Can deal with (start, stop] data and strata, as well as
#' sparse design matrices.
#'
#' WARNING: Users should not call \code{cox.fit} directly. Higher-level
#' functions in this package call \code{cox.fit} as a subroutine. If a
#' warm start object is provided, some of the other arguments in the function
#' may be overriden.
#'
#' \code{cox.fit} solves the elastic net problem for a single, user-specified
#' value of lambda. \code{cox.fit} works for Cox regression models, including
#' (start, stop] data and strata. It solves the problem using iteratively
#' reweighted least squares (IRLS). For each IRLS iteration, \code{cox.fit}
#' makes a quadratic (Newton) approximation of the log-likelihood, then calls
#' \code{elnet.fit} to minimize the resulting approximation.
#'
#' In terms of standardization: \code{cox.fit} does not standardize \code{x}
#' and \code{weights}. \code{penalty.factor} is standardized so that they sum
#' up to \code{nvars}.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed that any standardization needed
#' has already been done.
#' @param y Survival response variable, must be a Surv or stratifySurv object.
#' @param weights Observation weights. \code{cox.fit} does NOT standardize
#' these weights.
#' @param lambda A single value for the \code{lambda} hyperparameter.
#' @param alpha See glmnet help file
#' @param offset See glmnet help file
#' @param thresh Convergence threshold for coordinate descent. Each inner
#' coordinate-descent loop continues until the maximum change in the objective
#' after any coefficient update is less than thresh times the null deviance.
#' Default value is \code{1e-10}.
#' @param maxit Maximum number of passes over the data; default is \code{10^5}.
#' (If a warm start object is provided, the number of passes the warm start object
#' performed is included.)
#' @param penalty.factor See glmnet help file
#' @param exclude See glmnet help file
#' @param lower.limits See glmnet help file
#' @param upper.limits See glmnet help file
#' @param warm Either a \code{glmnetfit} object or a list (with name \code{beta}
#' containing coefficients) which can be used as a warm start. Default is
#' \code{NULL}, indicating no warm start. For internal use only.
#' @param from.cox.path Was \code{cox.fit()} called from \code{cox.path()}?
#' Default is FALSE.This has implications for computation of the penalty factors.
#' @param save.fit Return the warm start object? Default is FALSE.
#' @param trace.it Controls how much information is printed to screen. If
#' \code{trace.it=2}, some information about the fitting procedure is printed to
#' the console as the model is being fitted. Default is \code{trace.it=0}
#' (no information printed). (\code{trace.it=1} not used for compatibility with
#' \code{glmnet.path}.)
#'
#' @return An object with class "coxnet", "glmnetfit" and "glmnet". The list
#' returned contains more keys than that of a "glmnet" object.
#' \item{a0}{Intercept value, \code{NULL} for "cox" family.}
#' \item{beta}{A \code{nvars x 1} matrix of coefficients, stored in sparse matrix
#' format.}
#' \item{df}{The number of nonzero coefficients.}
#' \item{dim}{Dimension of coefficient matrix.}
#' \item{lambda}{Lambda value used.}
#' \item{dev.ratio}{The fraction of (null) deviance explained. The deviance
#' calculations incorporate weights if present in the model. The deviance is
#' defined to be 2*(loglike_sat - loglike), where loglike_sat is the log-likelihood
#' for the saturated model (a model with a free parameter per observation).
#' Hence dev.ratio=1-dev/nulldev.}
#' \item{nulldev}{Null deviance (per observation). This is defined to be
#' 2*(loglike_sat -loglike(Null)). The null model refers to the 0 model.}
#' \item{npasses}{Total passes over the data.}
#' \item{jerr}{Error flag, for warnings and errors (largely for internal
#' debugging).}
#' \item{offset}{A logical variable indicating whether an offset was included
#' in the model.}
#' \item{call}{The call that produced this object.}
#' \item{nobs}{Number of observations.}
#' \item{warm_fit}{If \code{save.fit=TRUE}, output of C++ routine, used for
#' warm starts. For internal use only.}
#' \item{family}{Family used for the model, always "cox".}
#' \item{converged}{A logical variable: was the algorithm judged to have
#' converged?}
#' \item{boundary}{A logical variable: is the fitted value on the boundary of
#' the attainable values?}
#' \item{obj_function}{Objective function value at the solution.}
cox.fit <- function(x, y, weights, lambda, alpha = 1.0, offset = rep(0, nobs),
                    thresh = 1e-10, maxit = 100000,
                    penalty.factor = rep(1.0, nvars), exclude = c(),
                    lower.limits = -Inf, upper.limits = Inf, warm = NULL,
                    from.cox.path = FALSE, save.fit = FALSE, trace.it = 0) {
  this.call <- match.call()
  control <- glmnet.control()

  ### Prepare all generic arguments
  nobs <- nrow(x)
  nvars <- ncol(x)
  is.offset <- !(missing(offset))
  if (is.offset == FALSE) {
    offset <- as.double(rep(0, nobs))
  }
  # add xm and xs attributes if they are missing for sparse x
  # glmnet.fit assumes that x is already standardized. Any standardization
  # the user wants should be done beforehand.
  if (inherits(x, "sparseMatrix")) {
    if ("xm" %in% names(attributes(x)) == FALSE)
      attr(x, "xm") <- rep(0.0, times = nvars)
    if ("xs" %in% names(attributes(x)) == FALSE)
      attr(x, "xs") <- rep(1.0, times = nvars)
  }

  # if calling from cox.path(), we do not need to check on exclude
  # and penalty.factor arguments as they have been prepared by cox.path()
  if (!from.cox.path) {
    # check and standardize penalty factors (to sum to nvars)
    if(any(penalty.factor == Inf)) {
      exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
      exclude = sort(unique(exclude))
    }
    if(length(exclude) > 0) {
      jd = match(exclude, seq(nvars), 0)
      if(!all(jd > 0)) stop ("Some excluded variables out of range")
      penalty.factor[jd] = 1 # ow can change lambda sequence
    }
    vp = pmax(0, penalty.factor)
    vp = as.double(vp * nvars / sum(vp))
  } else {
    vp <- as.double(penalty.factor)
  }

  ### check on limits
  lower.limits[lower.limits == -Inf] = -control$big
  upper.limits[upper.limits == Inf] = control$big
  if (length(lower.limits) < nvars)
    lower.limits = rep(lower.limits, nvars) else
      lower.limits = lower.limits[seq(nvars)]
  if (length(upper.limits) < nvars)
    upper.limits = rep(upper.limits, nvars) else
      upper.limits = upper.limits[seq(nvars)]
  ### end check on limits
  ### end preparation of generic arguments

  # compute null deviance
  if (is.null(warm)) {
    nulldev <- coxnet.deviance(y = y, offset = offset, weights = weights,
                               std.weights = FALSE)
    fit <- NULL
    coefold <- rep(0, nvars)   # initial coefs = 0
    eta <- offset
  } else {
    if (inherits(warm,"glmnetfit")) {
      if (!is(warm$warm_fit,"warmfit")) stop("Invalid warm start object")
      fit <- warm
      nulldev <- fit$nulldev
      coefold <- fit$warm_fit$a  # prev value for coefficients
      eta <- get_eta(x, coefold, 0) + offset
    } else if (inherits(warm,"list") && "beta" %in% names(warm)) {
      fit <- warm
      nulldev <- coxnet.deviance(y = y, offset = offset, weights = weights,
                                 std.weights = FALSE)
      coefold <- fit$beta   # prev value for coefficients
      eta <- get_eta(x, coefold, 0) + offset
      fit$a0 <- 0  # needed for compatibility with elnet.fit()
    } else {
      stop("Invalid warm start object")
    }
  }

  start <- NULL     # current value for coefficients
  obj_val_old <- cox_obj_function(y, eta, weights, lambda, alpha, coefold, vp)
  if (trace.it == 2) {
    cat("Warm Start Objective:", obj_val_old, fill = TRUE)
  }
  conv <- FALSE      # converged?

  # IRLS loop
  for (iter in 1L:control$mxitnr) {
    # compute working response and weights
    coxgrad_results <- coxgrad(eta, y, weights, std.weights = FALSE,
                               diag.hessian = TRUE)
    w <- -attributes(coxgrad_results)$diag_hessian
    z <- (eta - offset) - ifelse(w != 0, -coxgrad_results / w, 0)

    # have to update the weighted residual in our fit object
    # (in theory g and iy should be updated too, but we actually recompute g
    # and iy anyway in wls.f)
    if (!is.null(fit)) {
      fit$warm_fit$r <- w * (z - eta + offset)
    }

    # do WLS with warmstart from previous iteration
    fit <- elnet.fit(x, z, w, lambda, alpha, intercept = FALSE,
                     thresh = thresh, maxit = maxit, penalty.factor = vp,
                     exclude = exclude, lower.limits = lower.limits,
                     upper.limits = upper.limits, warm = fit,
                     from.glmnet.fit = TRUE, save.fit = TRUE)
    if (fit$jerr != 0) return(list(jerr = fit$jerr))

    # update coefficients, eta, mu and obj_val
    start <- fit$warm_fit$a
    eta <- get_eta(x, start, 0) + offset
    obj_val <- cox_obj_function(y, eta, weights, lambda, alpha, start, vp)
    if (trace.it == 2) cat("Iteration", iter, "Objective:", obj_val, fill = TRUE)

    boundary <- FALSE
    halved <- FALSE  # did we have to halve the step size?
    # if objective function is not finite, keep halving the stepsize until it is finite
    # for the halving step, we probably have to adjust fit$g as well?
    if (!is.finite(obj_val) || obj_val > control$big) {
      warning("Infinite objective function!", call. = FALSE)
      if (is.null(coefold))
        stop("no valid set of coefficients has been found: please supply starting values",
             call. = FALSE)
      warning("step size truncated due to divergence", call. = FALSE)
      ii <- 1
      while (!is.finite(obj_val) || obj_val > control$big) {
        if (ii > control$mxitnr)
          stop("inner loop 1; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold)/2
        eta <- get_eta(x, start, 0) + offset
        obj_val <- cox_obj_function(y, eta, weights, lambda, alpha, start, vp)
        if (trace.it == 2) cat("Iteration", iter, " Halved step 1, Objective:",
                               obj_val, fill = TRUE)
      }
      boundary <- TRUE
      halved <- TRUE
    }

    # if we did any halving, we have to update the coefficients, intercept
    # and weighted residual in the warm_fit object
    if (halved) {
      fit$warm_fit$a <- start
      fit$warm_fit$r <- w * (z - eta)
    }

    # test for convergence
    if (abs(obj_val - obj_val_old)/(0.1 + abs(obj_val)) < control$epsnr) {
      conv <- TRUE
      break
    }
    else {
      coefold <- start
      obj_val_old <- obj_val
    }
  }
  # end of IRLS loop

  # checks on convergence and fitted values
  if (!conv)
    warning("cox.fit: algorithm did not converge", call. = FALSE)

  # prepare output object
  if (save.fit == FALSE) {
    fit$warm_fit <- NULL
  }
  # overwrite values from elnet.fit object
  fit$a0 <- list(NULL)
  fit$call <- this.call
  fit$offset <- is.offset
  fit$nulldev <- nulldev
  fit$dev.ratio <- 1 - coxnet.deviance(y = y, pred = eta, weights = weights,
                                       std.weights = FALSE) / nulldev

  # add new key-value pairs to list
  fit$family <- "cox"
  fit$converged <- conv
  fit$boundary <- boundary
  fit$obj_function <- obj_val

  class(fit) <- c("coxnet", "glmnetfit", "glmnet")
  fit
}

#' Elastic net objective function value for Cox regression model
#'
#' Returns the elastic net objective function value for Cox regression model.
#'
#' @param y Survival response variable, must be a \code{Surv} or
#' \code{stratifySurv} object.
#' @param pred Model's predictions for \code{y}.
#' @param weights Observation weights.
#' @param lambda A single value for the \code{lambda} hyperparameter.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' @param coefficients The model's coefficients.
#' @param vp Penalty factors for each of the coefficients.
cox_obj_function <- function(y, pred, weights, lambda, alpha,
                             coefficients, vp) {
  coxnet.deviance(y = y, pred = pred, weights = weights, std.weights = FALSE) +
    lambda * pen_function(coefficients, alpha, vp)
}

#' Get lambda max for Cox regression model
#'
#' Return the lambda max value for Cox regression model, used for computing
#' initial lambda values. For internal use only.
#'
#' This function is called by \code{cox.path} for the value of lambda max.
#'
#' When \code{x} is not sparse, it is expected to already by centered and scaled.
#' When \code{x} is sparse, the function will get its attributes \code{xm} and
#' \code{xs} for its centering and scaling factors. The value of
#' \code{lambda_max} changes depending on whether \code{x} is centered and
#' scaled or not, so we need \code{xm} and \code{xs} to get the correct value.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed to be standardized.
#' @param y Survival response variable, must be a \code{Surv} or
#' \code{stratifySurv} object.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' @param weights Observation weights.
#' @param offset Offset for the model. Default is a zero vector of length
#' \code{nrow(y)}.
#' @param exclude Indices of variables to be excluded from the model.
#' @param vp Separate penalty factors can be applied to each coefficient.
get_cox_lambda_max <- function(x, y, alpha, weights = rep(1, nrow(x)),
                               offset = rep(0, nrow(x)), exclude = c(),
                               vp = rep(1, ncol(x))) {
  nobs <- nrow(x); nvars <- ncol(x)

  # extract strata (if any)
  if ("strata" %in% names(attributes(y))) {
    strata <- attr(y, "strata")
  } else {
    strata <- rep(1, nobs)
  }
  if (length(strata) != nobs) stop("length of strata != nobs")

  # if some penalty factors are zero, we need to compute eta
  vp_zero <- setdiff(which(vp == 0), exclude)
  if (length(vp_zero) > 0) {
    tempx <- x[, vp_zero, drop = FALSE]
    if(inherits(tempx, "sparseMatrix")) {
        attr(tempx, "xm") <- rep(0.0, length(vp_zero))
        attr(tempx, "xs") <- attr(x,"xs")[vp_zero]
      ## coxph cannot handle sparse x. Strata not needed because y is expected to be stratified
      fit <- cox.fit(x = tempx, y = y, offset = offset, weights = weights/sum(weights), lambda = 0)
      fit$beta <- fit$beta/attr(tempx,"xs")# need to put beta on the correct scale for next line to work
      eta <- as.numeric(predict(fit, newx = tempx, newoffset = offset, newstrata = strata))
    } else {
        eps <- glmnet.control()$epsnr

      if (length(unique(strata)) == 1) {
        fit <- survival::coxph(y ~ offset(offset) + tempx, weights = weights, eps = eps)
      } else {
        fit <- survival::coxph(y ~ offset(offset) + tempx + strata(strata),
                               weights = weights, eps = eps)
      }
      eta <- predict(fit, reference="sample") ## Coxph can do strata-specific centering
    }
  } else {
    eta <- offset
  }
  eta <- eta - mean(eta) ## keep numbers small; partial likelihood independent of centering
  ju <- rep(1, nvars)
  ju[exclude] <- 0 # we have already included constant variables in exclude

  # get cox gradient at "null" point
  # note that coxgrad already includes weights, so no need to include them
  # in subsequent computations
  null_grad <- coxgrad(eta, y, weights)

    if (inherits(x, "sparseMatrix")) {
        xm <- attr(x, "xm")
        xs <- attr(x, "xs")
        g <- abs((drop(t(null_grad) %*% x) - sum(null_grad) * xm) / xs)
    } else {
        g <- abs(drop(t(null_grad) %*% x))
    }
  g <- g / ifelse(vp > 0, vp, 1)
  g[ju == 0] <- 0
  lambda_max <- max(g) / max(alpha, 1e-3)
  return(lambda_max)
}

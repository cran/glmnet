#' Fit a GLM with elastic net regularization for a path of lambda values
#'
#' Fit a generalized linear model via penalized maximum likelihood for a path of
#' lambda values. Can deal with any GLM family.
#'
#' \code{glmnet.path} solves the elastic net problem for a path of lambda values.
#' It generalizes \code{glmnet::glmnet} in that it works for any GLM family.
#'
#' Sometimes the sequence is truncated before \code{nlambda} values of lambda
#' have been used. This happens when \code{glmnet.path} detects that the decrease
#' in deviance is marginal (i.e. we are near a saturated fit).
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. Can be a sparse matrix.
#' @param y Quantitative response variable.
#' @param weights Observation weights. Default is 1 for each observation.
#' @param lambda A user supplied lambda sequence. Typical usage is to have the
#' program compute its own lambda sequence based on \code{nlambda} and
#' \code{lambda.min.ratio}. Supplying a value of lambda overrides this.
#' @param nlambda The number of lambda values, default is 100.
#' @param lambda.min.ratio Smallest value for lambda as a fraction of lambda.max,
#' the (data derived) entry value (i.e. the smallest value for which all
#' coefficients are zero). The default depends on the sample size \code{nobs}
#' relative to the number of variables \code{nvars}. If \code{nobs >= nvars}, the
#' default is 0.0001, close to zero. If \code{nobs < nvars}, the default is 0.01.
#' A very small value of \code{lambda.min.ratio} will lead to a saturated fit
#' in the \code{nobs < nvars} case. This is undefined for some families of
#' models, and the function will exit gracefully when the percentage deviance
#' explained is almost 1.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.}
#' \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty.
#' @param offset A vector of length \code{nobs} that is included in the linear
#' predictor. Useful for the "poisson" family (e.g. log of exposure time), or
#' for refining a model by starting at a current fit. Default is NULL. If
#' supplied, then values must also be supplied to the \code{predict} function.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function. Default
#' is \code{gaussian()}. (See \code{\link[stats:family]{family}} for details on
#' family functions.)
#' @param standardize Logical flag for x variable standardization, prior to
#' fitting the model sequence. The coefficients are always returned on the
#' original scale. Default is \code{standardize=TRUE}. If variables are in the
#' same units already, you might not wish to standardize.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE)?
#' @param thresh Convergence threshold for coordinate descent. Each inner
#' coordinate-descent loop continues until the maximum change in the objective
#' after any coefficient update is less than thresh times the null deviance.
#' Default value is \code{1e-10}.
#' @param maxit Maximum number of passes over the data; default is \code{10^5}.
#' @param penalty.factor Separate penalty factors can be applied to each
#' coefficient. This is a number that multiplies \code{lambda} to allow differential
#' shrinkage. Can be 0 for some variables, which implies no shrinkage, and that
#' variable is always included in the model. Default is 1 for all variables (and
#' implicitly infinity for variables listed in exclude). Note: the penalty
#' factors are internally rescaled to sum to \code{nvars}.
#' @param exclude Indices of variables to be excluded from the model. Default is
#' none. Equivalent to an infinite penalty factor.
#' @param lower.limits Vector of lower limits for each coefficient; default
#' \code{-Inf}. Each of these must be non-positive. Can be presented as a single
#' value (which will then be replicated), else a vector of length \code{nvars}.
#' @param upper.limits Vector of upper limits for each coefficient; default
#' \code{Inf}. See \code{lower.limits}.
#' @param trace.it Controls how much information is printed to screen. Default is
#' \code{trace.it=0} (no information printed). If \code{trace.it=1}, a progress
#' bar is displayed. If \code{trace.it=2}, some information about the fitting
#' procedure is printed to the console as the model is being fitted.
#'
#' @return An object with class "glmnetfit" and "glmnet".
#' \item{a0}{Intercept sequence of length \code{length(lambda)}.}
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
#' 2*(loglike_sat -loglike(Null)). The null model refers to the intercept model.}
#' \item{npasses}{Total passes over the data summed over all lambda values.}
#' \item{jerr}{Error flag, for warnings and errors (largely for internal
#' debugging).}
#' \item{offset}{A logical variable indicating whether an offset was included
#' in the model.}
#' \item{call}{The call that produced this object.}
#' \item{family}{Family used for the model.}
#' \item{nobs}{Number of observations.}
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100 * 20), nrow = 100)
#' y <- ifelse(rnorm(100) > 0, 1, 0)
#'
#' # binomial with probit link
#' fit1 <- glmnet:::glmnet.path(x, y, family = binomial(link = "probit"))
glmnet.path <- function(x, y, weights=NULL, lambda = NULL, nlambda = 100,
                        lambda.min.ratio = ifelse(nobs<nvars, 0.01, 0.0001),
                        alpha = 1.0, offset = NULL, family = gaussian(),
                        standardize = TRUE, intercept = TRUE, thresh = 1e-10, maxit = 100000,
                        penalty.factor = rep(1.0, nvars), exclude = integer(0), lower.limits = -Inf,
                        upper.limits = Inf, trace.it = 0) {

    ### Check on family argument
    if(is.function(family))family=family()
    if(!inherits(family,"family"))
        stop("Invalid family argument; must be either character, function or family object")

    ### Prepare all the generic arguments
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
    if(is.null(np) || (np[2] <= 1)) stop("x should be a matrix with 2 or more columns")
    nobs = as.integer(np[1]); nvars = as.integer(np[2])

    # get feature variable names
    vnames=colnames(x)
    if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")

    # check weights
    if(is.null(weights)) weights = rep(1,nobs)
    else if (length(weights) != nobs)
        stop(paste("Number of elements in weights (",length(weights),
                   ") not equal to the number of rows of x (",nobs,")",sep=""))
    weights <- as.double(weights)

    ## initialize from family function. Makes y a vector in case of binomial, and possibly changes weights
    ## Expects nobs to be defined, and creates n and mustart (neither used here)
    ## Some cases expect to see things, so we set it up just to make it work
    etastart=0;mustart=NULL;start=NULL
    eval(family$initialize)
    ##
    ## Just in case this was not done in initialize()
    y <- drop(y)  # we don't like matrix responses

    is.offset <- !(is.null(offset))
    if (is.offset == FALSE) {
        offset <- as.double(y * 0) #keeps the shape of y
    }
    # infinite penalty factor vars are excluded
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
    if (intercept) {
        xm <- meansd$mean
    } else {
        xm <- rep(0.0, times = nvars)
    }
    if (standardize) {
        xs <- meansd$sd
    } else {
        xs <- rep(1.0, times = nvars)
    }
    if (!inherits(x, "sparseMatrix")) {
        x <- scale(x, xm, xs)
    } else {
        attr(x, "xm") <- xm
        attr(x, "xs") <- xs
    }
    lower.limits <- lower.limits * xs
    upper.limits <- upper.limits * xs

    # get null deviance and lambda max
    start_val <- get_start(x, y, weights, family, intercept, is.offset,
                           offset, exclude, vp, alpha)

    # work out lambda values
    nlam = as.integer(nlambda)
    user_lambda = FALSE   # did user provide their own lambda values?
    if (is.null(lambda)) {
        if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")

        # compute lambda max: to add code here
        lambda_max <- start_val$lambda_max

        # compute lambda sequence
        ulam <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                        length.out = nlam))
    } else {  # user provided lambda values
        user_lambda = TRUE
        if (any(lambda < 0)) stop("lambdas should be non-negative")
        ulam = as.double(rev(sort(lambda)))
        nlam = as.integer(length(lambda))
    }

    # start progress bar
    if (trace.it == 1) pb <- utils::txtProgressBar(min = 0, max = nlam, style = 3)

    a0 <- rep(NA, length = nlam)
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
        fit <- glmnet.fit(x, y, weights / sum(weights), cur_lambda, alpha = alpha, offset = offset,
                          family = family, intercept = intercept, thresh = thresh,
                          maxit = maxit, penalty.factor = vp, exclude = exclude,
                          lower.limits = lower.limits, upper.limits = upper.limits,
                          warm = fit, from.glmnet.path = TRUE, save.fit = TRUE,
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

        a0[k] <- fit$a0
        beta[, k] <- as.matrix(fit$beta)
        dev.ratio[k] <- fit$dev.ratio

        # early stopping if dev.ratio almost 1 or no improvement
        if (k >= mnl && user_lambda == FALSE) {
            if (dev.ratio[k] > control$devmax) break
            else if (k > 1) {
                if (family$family == "gaussian") {
                    if (dev.ratio[k] - dev.ratio[k-1] < control$fdev * dev.ratio[k])
                        break
                } else if (family$family == "poisson") {
                    if (dev.ratio[k] - dev.ratio[k - mnl + 1] <
                        10 * control$fdev * dev.ratio[k])
                        break
                } else if (dev.ratio[k] - dev.ratio[k-1] < control$fdev) break
            }
        }
    }
    if (trace.it == 1) {
        utils::setTxtProgressBar(pb, nlam)
        cat("", fill = TRUE)
    }

    # truncate a0, beta, dev.ratio, lambda if necessary
    if (k < nlam) {
        a0 <- a0[1:k]
        beta <- beta[, 1:k, drop = FALSE]
        dev.ratio <- dev.ratio[1:k]
        ulam <- ulam[1:k]
    }

    # return coefficients to original scale (because of x standardization)
    beta <- beta / xs
    a0 <- a0 - colSums(beta * xm)

    # output
    stepnames <- paste0("s", 0:(length(ulam) - 1))
    out <- list()
    out$a0 <- a0
    names(out$a0) <- stepnames
    out$beta <- Matrix::Matrix(beta, sparse = TRUE,
                               dimnames = list(vnames, stepnames))
    out$df <- colSums(abs(beta) > 0)
    out$dim <- dim(beta)
    out$lambda <- ulam
    out$dev.ratio <- dev.ratio
    out$nulldev <- start_val$nulldev
    out$npasses <- fit$npasses
    out$jerr <- fit$jerr
    out$offset <- is.offset
    out$call <- this.call
    out$family <- family
    out$nobs <- nobs
    class(out) <- c("glmnetfit", "glmnet")

    return(out)
}
#' Fit a GLM with elastic net regularization for a single value of lambda
#'
#' Fit a generalized linear model via penalized maximum likelihood for a single
#' value of lambda. Can deal with any GLM family.
#'
#' WARNING: Users should not call \code{glmnet.fit} directly. Higher-level functions
#' in this package call \code{glmnet.fit} as a subroutine. If a warm start object
#' is provided, some of the other arguments in the function may be overriden.
#'
#' \code{glmnet.fit} solves the elastic net problem for a single, user-specified
#' value of lambda. \code{glmnet.fit} works for any GLM family. It solves the
#' problem using iteratively reweighted least squares (IRLS). For each IRLS
#' iteration, \code{glmnet.fit} makes a quadratic (Newton) approximation of the
#' log-likelihood, then calls \code{elnet.fit} to minimize the resulting
#' approximation.
#'
#' In terms of standardization: \code{glmnet.fit} does not standardize \code{x}
#' and \code{weights}. \code{penalty.factor} is standardized so that they sum up
#' to \code{nvars}.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed that any standardization needed
#' has already been done.
#' @param y Quantitative response variable.
#' @param weights Observation weights. \code{glmnet.fit} does NOT standardize
#' these weights.
#' @param lambda A single value for the \code{lambda} hyperparameter.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.}
#' \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty.
#' @param offset A vector of length \code{nobs} that is included in the linear
#' predictor. Useful for the "poisson" family (e.g. log of exposure time), or
#' for refining a model by starting at a current fit. Default is NULL. If
#' supplied, then values must also be supplied to the \code{predict} function.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function. Default
#' is \code{gaussian()}. (See \code{\link[stats:family]{family}} for details on
#' family functions.)
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE)?
#' @param thresh Convergence threshold for coordinate descent. Each inner
#' coordinate-descent loop continues until the maximum change in the objective
#' after any coefficient update is less than thresh times the null deviance.
#' Default value is \code{1e-10}.
#' @param maxit Maximum number of passes over the data; default is \code{10^5}.
#' (If a warm start object is provided, the number of passes the warm start object
#' performed is included.)
#' @param penalty.factor Separate penalty factors can be applied to each
#' coefficient. This is a number that multiplies \code{lambda} to allow differential
#' shrinkage. Can be 0 for some variables, which implies no shrinkage, and that
#' variable is always included in the model. Default is 1 for all variables (and
#' implicitly infinity for variables listed in exclude). Note: the penalty
#' factors are internally rescaled to sum to \code{nvars}.
#' @param exclude Indices of variables to be excluded from the model. Default is
#' none. Equivalent to an infinite penalty factor.
#' @param lower.limits Vector of lower limits for each coefficient; default
#' \code{-Inf}. Each of these must be non-positive. Can be presented as a single
#' value (which will then be replicated), else a vector of length \code{nvars}.
#' @param upper.limits Vector of upper limits for each coefficient; default
#' \code{Inf}. See \code{lower.limits}.
#' @param warm Either a \code{glmnetfit} object or a list (with names \code{beta}
#' and \code{a0} containing coefficients and intercept respectively) which can
#' be used as a warm start. Default is \code{NULL}, indicating no warm start.
#' For internal use only.
#' @param from.glmnet.path Was \code{glmnet.fit()} called from \code{glmnet.path()}?
#' Default is FALSE.This has implications for computation of the penalty factors.
#' @param save.fit Return the warm start object? Default is FALSE.
#' @param trace.it Controls how much information is printed to screen. If
#' \code{trace.it=2}, some information about the fitting procedure is printed to
#' the console as the model is being fitted. Default is \code{trace.it=0}
#' (no information printed). (\code{trace.it=1} not used for compatibility with
#' \code{glmnet.path}.)
#'
#' @return An object with class "glmnetfit" and "glmnet". The list
#' returned contains more keys than that of a "glmnet" object.
#' \item{a0}{Intercept value.}
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
#' 2*(loglike_sat -loglike(Null)). The null model refers to the intercept model.}
#' \item{npasses}{Total passes over the data.}
#' \item{jerr}{Error flag, for warnings and errors (largely for internal
#' debugging).}
#' \item{offset}{A logical variable indicating whether an offset was included
#' in the model.}
#' \item{call}{The call that produced this object.}
#' \item{nobs}{Number of observations.}
#' \item{warm_fit}{If \code{save.fit=TRUE}, output of C++ routine, used for
#' warm starts. For internal use only.}
#' \item{family}{Family used for the model.}
#' \item{converged}{A logical variable: was the algorithm judged to have
#' converged?}
#' \item{boundary}{A logical variable: is the fitted value on the boundary of
#' the attainable values?}
#' \item{obj_function}{Objective function value at the solution.}
#'
glmnet.fit <- function(x, y, weights, lambda, alpha = 1.0,
                       offset = rep(0, nobs), family = gaussian(),
                       intercept = TRUE, thresh = 1e-10, maxit = 100000,
                       penalty.factor = rep(1.0, nvars), exclude = c(), lower.limits = -Inf,
                       upper.limits = Inf, warm = NULL, from.glmnet.path = FALSE,
                       save.fit = FALSE, trace.it = 0) {
    this.call <- match.call()
    control <- glmnet.control()

    ### Prepare all the generic arguments
    nobs <- nrow(x)
    nvars <- ncol(x)
    is.offset <- !(missing(offset))
    if (is.offset == FALSE) {
        offset <- as.double(y * 0) #keeps the shape of y
    } else if (is.null(offset)) {
        offset <- rep(0, nobs)
        is.offset = FALSE
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

    # if calling from glmnet.path(), we do not need to check on exclude
    # and penalty.factor arguments as they have been prepared by glmnet.path()
    if (!from.glmnet.path) {
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

    # get the relevant family functions
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
             call. = FALSE)
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x))
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)

    # computation of null deviance (get mu in the process)
    if (is.null(warm)) {
        start_val <- get_start(x, y, weights, family, intercept, is.offset,
                               offset, exclude, vp, alpha)
        nulldev <- start_val$nulldev
        mu <- start_val$mu
        fit <- NULL
        coefold <- rep(0, nvars)   # initial coefs = 0
        eta <- family$linkfun(mu)
        intold <- (eta - offset)[1]
    } else {
        if (inherits(warm,"glmnetfit")) {
            if (!is(warm$warm_fit,"warmfit")) stop("Invalid warm start object")
            fit <- warm
            nulldev <- fit$nulldev
            coefold <- fit$warm_fit$a   # prev value for coefficients
            intold <- fit$warm_fit$aint    # prev value for intercept
            eta <- get_eta(x, coefold, intold)
            mu <- linkinv(eta <- eta + offset)
        } else if (inherits(warm,"list") && "a0" %in% names(warm) &&
                   "beta" %in% names(warm)) {
            nulldev <- get_start(x, y, weights, family, intercept, is.offset,
                                   offset, exclude, vp, alpha)$nulldev
            fit <- warm
            coefold <- fit$beta   # prev value for coefficients
            intold <- fit$a0    # prev value for intercept
            eta <- get_eta(x, coefold, intold)
            mu <- linkinv(eta <- eta + offset)
        } else {
            stop("Invalid warm start object")
        }
    }

    if (!(validmu(mu) && valideta(eta)))
        stop("cannot find valid starting values: please specify some",
             call. = FALSE)

    start <- NULL     # current value for coefficients
    start_int <- NULL # current value for intercept
    obj_val_old <- obj_function(y, mu, weights, family, lambda, alpha, coefold, vp)
    if (trace.it == 2) {
        cat("Warm Start Objective:", obj_val_old, fill = TRUE)
    }
    conv <- FALSE      # converged?

    # IRLS loop
    for (iter in 1L:control$mxitnr) {
        # some checks for NAs/zeros
        varmu <- variance(mu)
        if (anyNA(varmu)) stop("NAs in V(mu)")
        if (any(varmu == 0)) stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val))) stop("NAs in d(mu)/d(eta)")

        # compute working response and weights
        z <- (eta - offset) + (y - mu)/mu.eta.val
        w <- (weights * mu.eta.val^2)/variance(mu)

        # have to update the weighted residual in our fit object
        # (in theory g and iy should be updated too, but we actually recompute g
        # and iy anyway in wls.f)
        if (!is.null(fit)) {
            fit$warm_fit$r <- w * (z - eta + offset)
        }

        # do WLS with warmstart from previous iteration
        fit <- elnet.fit(x, z, w, lambda, alpha, intercept,
                         thresh = thresh, maxit = maxit, penalty.factor = vp,
                         exclude = exclude, lower.limits = lower.limits,
                         upper.limits = upper.limits, warm = fit,
                         from.glmnet.fit = TRUE, save.fit = TRUE)
        if (fit$jerr != 0) return(list(jerr = fit$jerr))

        # update coefficients, eta, mu and obj_val
        start <- fit$warm_fit$a
        start_int <- fit$warm_fit$aint
        eta <- get_eta(x, start, start_int)
        mu <- linkinv(eta <- eta + offset)
        obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp)
        if (trace.it == 2) cat("Iteration", iter, "Objective:", obj_val, fill = TRUE)

        boundary <- FALSE
        halved <- FALSE  # did we have to halve the step size?
        # if objective function is not finite, keep halving the stepsize until it is finite
        # for the halving step, we probably have to adjust fit$g as well?
        if (!is.finite(obj_val) || obj_val > control$big) {
            warning("Infinite objective function!", call. = FALSE)
            if (is.null(coefold) || is.null(intold))
                stop("no valid set of coefficients has been found: please supply starting values",
                     call. = FALSE)
            warning("step size truncated due to divergence", call. = FALSE)
            ii <- 1
            while (!is.finite(obj_val) || obj_val > control$big) {
                if (ii > control$mxitnr)
                    stop("inner loop 1; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold)/2
                start_int <- (start_int + intold)/2
                eta <- get_eta(x, start, start_int)
                mu <- linkinv(eta <- eta + offset)
                obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp)
                if (trace.it == 2) cat("Iteration", iter, " Halved step 1, Objective:",
                               obj_val, fill = TRUE)
            }
            boundary <- TRUE
            halved <- TRUE
        }
        # if some of the new eta or mu are invalid, keep halving stepsize until valid
        if (!(valideta(eta) && validmu(mu))) {
            warning("Invalid eta/mu!", call. = FALSE)
            if (is.null(coefold) || is.null(intold))
                stop("no valid set of coefficients has been found: please supply starting values",
                     call. = FALSE)
            warning("step size truncated: out of bounds", call. = FALSE)
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$mxitnr)
                    stop("inner loop 2; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold)/2
                start_int <- (start_int + intold)/2
                eta <- get_eta(x, start, start_int)
                mu <- linkinv(eta <- eta + offset)
            }
            boundary <- TRUE
            halved <- TRUE
            obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp)
            if (trace.it == 2) cat("Iteration", iter, " Halved step 2, Objective:", obj_val, fill = TRUE)
        }
        # extra halving step if objective function value actually increased
        if (obj_val > obj_val_old + 1e-7) {
            ii <- 1
            while (obj_val > obj_val_old + 1e-7) {
                if (ii > control$mxitnr)
                    stop("inner loop 3; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold)/2
                start_int <- (start_int + intold)/2
                eta <- get_eta(x, start, start_int)
                mu <- linkinv(eta <- eta + offset)
                obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp)
                if (trace.it == 2) cat("Iteration", iter, " Halved step 3, Objective:",
                                       obj_val, fill = TRUE)
            }
            halved <- TRUE
        }

        # if we did any halving, we have to update the coefficients, intercept
        # and weighted residual in the warm_fit object
        if (halved) {
            fit$warm_fit$a <- start
            fit$warm_fit$aint <- start_int
            fit$warm_fit$r <- w * (z - eta)
        }

        # test for convergence
        if (abs(obj_val - obj_val_old)/(0.1 + abs(obj_val)) < control$epsnr) {
            conv <- TRUE
            break
        }
        else {
            coefold <- start
            intold <- start_int
            obj_val_old <- obj_val
        }
    }
    # end of IRLS loop

    # checks on convergence and fitted values
    if (!conv)
        warning("glmnet.fit: algorithm did not converge", call. = FALSE)
    if (boundary)
        warning("glmnet.fit: algorithm stopped at boundary value", call. = FALSE)

    # some extra warnings, printed only if trace.it == 2
    if (trace.it == 2) {
        eps <- 10 * .Machine$double.eps
        if ((family$family == "binomial") && (any(mu > 1 - eps) || any(mu < eps)))
                warning("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                        call. = FALSE)
        if ((family$family == "poisson") && (any(mu < eps)))
                warning("glm.fit: fitted rates numerically 0 occurred",
                        call. = FALSE)
    }

    # prepare output object
    if (save.fit == FALSE) {
        fit$warm_fit <- NULL
    }
    # overwrite values from elnet.fit object
    fit$call <- this.call
    fit$offset <- is.offset
    fit$nulldev <- nulldev
    fit$dev.ratio <- 1 - dev_function(y, mu, weights, family) / nulldev

    # add new key-value pairs to list
    fit$family <- family
    fit$converged <- conv
    fit$boundary <- boundary
    fit$obj_function <- obj_val

    class(fit) <- c("glmnetfit", "glmnet")
    fit
}

#' Solve weighted least squares (WLS) problem for a single lambda value
#'
#' Solves the weighted least squares (WLS) problem for a single lambda value. Internal
#' function that users should not call directly.
#'
#' WARNING: Users should not call \code{elnet.fit} directly. Higher-level functions
#' in this package call \code{elnet.fit} as a subroutine. If a warm start object
#' is provided, some of the other arguments in the function may be overriden.
#'
#' \code{elnet.fit} is essentially a wrapper around a C++ subroutine which
#' minimizes
#'
#' \deqn{1/2 \sum w_i (y_i - X_i^T \beta)^2 + \sum \lambda \gamma_j
#' [(1-\alpha)/2 \beta^2+\alpha|\beta|],}
#'
#' over \eqn{\beta}, where \eqn{\gamma_j} is the relative penalty factor on the
#' jth variable. If \code{intercept = TRUE}, then the term in the first sum is
#' \eqn{w_i (y_i - \beta_0 - X_i^T \beta)^2}, and we are minimizing over both
#' \eqn{\beta_0} and \eqn{\beta}.
#'
#' None of the inputs are standardized except for \code{penalty.factor}, which
#' is standardized so that they sum up to \code{nvars}.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed that any standardization needed
#' has already been done.
#' @param y Quantitative response variable.
#' @param weights Observation weights. \code{elnet.fit} does NOT standardize
#' these weights.
#' @param lambda A single value for the \code{lambda} hyperparameter.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.}
#' \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE)?
#' @param thresh Convergence threshold for coordinate descent. Each inner
#' coordinate-descent loop continues until the maximum change in the objective
#' after any coefficient update is less than thresh times the null deviance.
#' Default value is \code{1e-7}.
#' @param maxit Maximum number of passes over the data; default is \code{10^5}.
#' (If a warm start object is provided, the number of passes the warm start object
#' performed is included.)
#' @param penalty.factor Separate penalty factors can be applied to each
#' coefficient. This is a number that multiplies \code{lambda} to allow differential
#' shrinkage. Can be 0 for some variables, which implies no shrinkage, and that
#' variable is always included in the model. Default is 1 for all variables (and
#' implicitly infinity for variables listed in exclude). Note: the penalty
#' factors are internally rescaled to sum to \code{nvars}.
#' @param exclude Indices of variables to be excluded from the model. Default is
#' none. Equivalent to an infinite penalty factor.
#' @param lower.limits Vector of lower limits for each coefficient; default
#' \code{-Inf}. Each of these must be non-positive. Can be presented as a single
#' value (which will then be replicated), else a vector of length \code{nvars}.
#' @param upper.limits Vector of upper limits for each coefficient; default
#' \code{Inf}. See \code{lower.limits}.
#' @param warm Either a \code{glmnetfit} object or a list (with names \code{beta}
#' and \code{a0} containing coefficients and intercept respectively) which can
#' be used as a warm start. Default is \code{NULL}, indicating no warm start.
#' For internal use only.
#' @param from.glmnet.fit Was \code{elnet.fit()} called from \code{glmnet.fit()}?
#' Default is FALSE.This has implications for computation of the penalty factors.
#' @param save.fit Return the warm start object? Default is FALSE.
#'
#' @return An object with class "glmnetfit" and "glmnet". The list returned has
#' the same keys as that of a \code{glmnet} object, except that it might have an
#' additional \code{warm_fit} key.
#' \item{a0}{Intercept value.}
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
#' 2*(loglike_sat -loglike(Null)). The null model refers to the intercept model.}
#' \item{npasses}{Total passes over the data.}
#' \item{jerr}{Error flag, for warnings and errors (largely for internal
#' debugging).}
#' \item{offset}{Always FALSE, since offsets do not appear in the WLS problem.
#' Included for compability with glmnet output.}
#' \item{call}{The call that produced this object.}
#' \item{nobs}{Number of observations.}
#' \item{warm_fit}{If \code{save.fit=TRUE}, output of C++ routine, used for
#' warm starts. For internal use only.}
#'
elnet.fit <- function(x, y, weights, lambda, alpha = 1.0, intercept = TRUE,
                      thresh = 1e-7, maxit = 100000,
                      penalty.factor = rep(1.0, nvars), exclude = c(),
                      lower.limits = -Inf, upper.limits = Inf, warm = NULL,
                      from.glmnet.fit = FALSE, save.fit = FALSE) {
    this.call <- match.call()
    internal.parms <- glmnet.control()

    # compute null deviance
    ybar <- weighted.mean(y, weights)
    nulldev <- sum(weights * (y - ybar)^2)

    # if class "glmnetfit" warmstart object provided, pull whatever we want out of it
    # else, prepare arguments, then check if coefs provided as warmstart
    # (if only coefs are given as warmstart, we prepare the other arguments
    # as if no warmstart was provided)
    if (!is.null(warm) && "glmnetfit" %in% class(warm)) {
        warm <- warm$warm_fit
        if (!is(warm,"warmfit")) stop("Invalid warm start object")

        a <- warm$a
        aint <- warm$aint
        alm0 <- warm$almc
        cl <- warm$cl
        g <- warm$g
        ia <- warm$ia
        iy <- warm$iy
        iz <- warm$iz
        ju <- warm$ju
        m <- warm$m
        mm <- warm$mm
        nino <- warm$nino
        nobs <- warm$no
        nvars <- warm$ni
        nlp <- warm$nlp
        nx <- warm$nx
        r <- warm$r
        rsqc <- warm$rsqc
        xv <- warm$xv
        vp <- warm$vp
    } else {
        nobs <- as.integer(nrow(x))
        nvars <- as.integer(ncol(x))

        # if calling from glmnet.fit(), we do not need to check on exclude
        # and penalty.factor arguments as they have been prepared by glmnet.fit()
        # Also exclude will include variance 0 columns
        if (!from.glmnet.fit) {
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
        # compute ju
        # assume that there are no constant variables
        ju <- rep(1, nvars)
        ju[exclude] <- 0
        ju <- as.integer(ju)

        # compute cl from lower.limits and upper.limits
        lower.limits[lower.limits == -Inf] <- -internal.parms$big
        upper.limits[upper.limits == Inf] <- internal.parms$big
        if (length(lower.limits) < nvars)
            lower.limits = rep(lower.limits, nvars) else
                lower.limits = lower.limits[seq(nvars)]
        if (length(upper.limits) < nvars)
            upper.limits = rep(upper.limits, nvars) else
                upper.limits = upper.limits[seq(nvars)]
        cl <- rbind(lower.limits, upper.limits)
        storage.mode(cl) = "double"

        nx <- as.integer(nvars)

        a <- double(nvars)
        aint <- double(1)
        alm0 <- double(1)
        g <- double(nvars)
        ia <- integer(nx)
        iy <- integer(nvars)
        iz <- integer(1)
        m <- as.integer(1)
        mm <- integer(nvars)
        nino <- integer(1)
        nlp <- integer(1)
        r <- weights * y
        rsqc <- double(1)
        xv <- double(nvars)

        # check if coefs were provided as warmstart: if so, use them
        if (!is.null(warm)) {
            if (inherits(warm,"list") && "a0" %in% names(warm) &&
                "beta" %in% names(warm)) {
                a <- as.double(warm$beta)
                aint <- as.double(warm$a0)
                mu <- drop(x %*% a + aint)
                r <- weights * (y - mu)
                rsqc <- 1 - sum(weights * (y - mu)^2) / nulldev
            } else {
                stop("Invalid warm start object")
            }
        }
    }

    # for the parameters here, we are overriding the values provided by the
    # warmstart object
    alpha <- as.double(alpha)
    almc <- as.double(lambda)
    intr <- as.integer(intercept)
    jerr <- integer(1)
    maxit <- as.integer(maxit)
    thr <- as.double(thresh)
    v <- as.double(weights)

    a.new <- a
    a.new[0] <- a.new[0] # induce a copy

    # take out components of x and run C++ subroutine
    if (inherits(x, "sparseMatrix")) {
        xm <- as.double(attr(x, "xm"))
        xs <- as.double(attr(x, "xs"))
        wls_fit <- spwls_exp(alm0=alm0,almc=almc,alpha=alpha,m=m,no=nobs,ni=nvars,
                             x=x,xm=xm,xs=xs,r=r,xv=xv,v=v,intr=intr,ju=ju,vp=vp,cl=cl,nx=nx,thr=thr,
                             maxit=maxit,a=a.new,aint=aint,g=g,ia=ia,iy=iy,iz=iz,mm=mm,
                             nino=nino,rsqc=rsqc,nlp=nlp,jerr=jerr)
    } else {
        wls_fit <- wls_exp(alm0=alm0,almc=almc,alpha=alpha,m=m,no=nobs,ni=nvars,
                           x=x,r=r,xv=xv,v=v,intr=intr,ju=ju,vp=vp,cl=cl,nx=nx,thr=thr,
                           maxit=maxit,a=a.new,aint=aint,g=g,ia=ia,iy=iy,iz=iz,mm=mm,
                           nino=nino,rsqc=rsqc,nlp=nlp,jerr=jerr)
    }

    # if error code > 0, fatal error occurred: stop immediately
    # if error code < 0, non-fatal error occurred: return error code
    if (wls_fit$jerr > 0) {
        errmsg <- jerr.glmnetfit(wls_fit$jerr, maxit)
        stop(errmsg$msg, call. = FALSE)
    } else if (wls_fit$jerr < 0)
        return(list(jerr = wls_fit$jerr))
    warm_fit <- wls_fit[c("almc", "r", "xv", "ju", "vp", "cl", "nx",
                          "a", "aint", "g", "ia", "iy", "iz", "mm", "nino",
                          "rsqc", "nlp")]
    warm_fit[['m']] <- m
    warm_fit[['no']] <- nobs
    warm_fit[['ni']] <- nvars
    class(warm_fit) <- "warmfit"

    beta <- Matrix::Matrix(wls_fit$a, sparse = TRUE)

    out <- list(a0 = wls_fit$aint, beta = beta, df = sum(abs(beta) > 0),
                dim = dim(beta), lambda = lambda, dev.ratio = wls_fit$rsqc,
                nulldev = nulldev, npasses = wls_fit$nlp, jerr = wls_fit$jerr,
                offset = FALSE, call = this.call, nobs = nobs, warm_fit = warm_fit)
    if (save.fit == FALSE) {
        out$warm_fit <- NULL
    }
    class(out) <- c("glmnetfit", "glmnet")
    out
}

#' Get null deviance, starting mu and lambda max
#'
#' Return the null deviance, starting mu and lambda max values for
#' initialization. For internal use only.
#'
#' This function is called by \code{glmnet.path} for null deviance, starting mu
#' and lambda max values. It is also called by \code{glmnet.fit} when used
#' without warmstart, but they only use the null deviance and starting mu values.
#'
#' When \code{x} is not sparse, it is expected to already by centered and scaled.
#' When \code{x} is sparse, the function will get its attributes \code{xm} and
#' \code{xs} for its centering and scaling factors.
#'
#' Note that whether \code{x} is centered & scaled or not, the values of \code{mu}
#' and \code{nulldev} don't change. However, the value of \code{lambda_max} does
#' change, and we need \code{xm} and \code{xs} to get the correct value.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed to be standardized.
#' @param y Quantitative response variable.
#' @param weights Observation weights.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function.
#' (See \code{\link[stats:family]{family}} for details on family functions.)
#' @param intercept Does the model we are fitting have an intercept term or not?
#' @param is.offset Is the model being fit with an offset or not?
#' @param offset Offset for the model. If \code{is.offset=FALSE}, this should be
#' a zero vector of the same length as \code{y}.
#' @param exclude Indices of variables to be excluded from the model.
#' @param vp Separate penalty factors can be applied to each coefficient.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
get_start <- function(x, y, weights, family, intercept, is.offset, offset,
                      exclude, vp, alpha) {
    nobs <- nrow(x); nvars <- ncol(x)

    # compute mu and null deviance
    # family = binomial() gives us warnings due to non-integer weights
    # to avoid, suppress warnings
    if (intercept) {
        if (is.offset) {
            suppressWarnings(tempfit <- glm(y ~ 1, family = family,
                                            weights = weights, offset = offset))
            mu <- tempfit$fitted.values
        } else {
            mu <- rep(weighted.mean(y,weights), times = nobs)
        }
    } else {
        mu <- family$linkinv(offset)
    }
    nulldev <- dev_function(y, mu, weights, family)

    # if some penalty factors are zero, we have to recompute mu
    vp_zero <- setdiff(which(vp == 0), exclude)
    if (length(vp_zero) > 0) {
        tempx <- x[, vp_zero, drop = FALSE]
        if (inherits(tempx, "sparseMatrix")) {
            tempfit <- glmnet.fit(tempx, y, intercept = intercept, family = family,
                                  weights = weights/sum(weights), offset = offset, lambda = 0)
            mu <- predict(tempfit, newx=tempx, newoffset=offset, type = "response")
        } else {
            if (intercept) {
                tempx <- cbind(1,tempx)
            }
            tempfit <- glm.fit(tempx, y, family = family, weights = weights, offset = offset)
            mu <- tempfit$fitted.values
        }
    }
    # compute lambda max
    ju <- rep(1, nvars)
    ju[exclude] <- 0 # we have already included constant variables in exclude
    r <- y - mu
    eta <- family$linkfun(mu)
    v <- family$variance(mu)
    m.e <- family$mu.eta(eta)
    weights <- weights / sum(weights)
    rv <- r / v * m.e * weights
    if (inherits(x, "sparseMatrix")) {
        xm <- attr(x, "xm")
        xs <- attr(x, "xs")
        g <- abs((drop(t(rv) %*% x) - sum(rv) * xm) / xs)
    } else {
        g <- abs(drop(t(rv) %*% x))
    }
    g <- g * ju / ifelse(vp > 0, vp, 1)
    lambda_max <- max(g) / max(alpha, 1e-3)

    list(nulldev = nulldev, mu = mu, lambda_max = lambda_max)
}
#' Elastic net objective function value
#'
#' Returns the elastic net objective function value.
#'
#' @param y Quantitative response variable.
#' @param mu Model's predictions for \code{y}.
#' @param weights Observation weights.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function.
#' @param lambda A single value for the \code{lambda} hyperparameter.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' @param coefficients The model's coefficients (excluding intercept).
#' @param vp Penalty factors for each of the coefficients.
obj_function <- function(y, mu, weights, family,
                         lambda, alpha, coefficients, vp) {
    dev_function(y, mu, weights, family) / 2 +
        lambda * pen_function(coefficients, alpha, vp)
}

#' Elastic net penalty value
#'
#' Returns the elastic net penalty value without the \code{lambda} factor.
#'
#' The penalty is defined as
#' \deqn{(1-\alpha)/2 \sum vp_j \beta_j^2 + \alpha \sum vp_j |\beta|.}
#' Note the omission of the multiplicative \code{lambda} factor.
#'
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' @param coefficients The model's coefficients (excluding intercept).
#' @param vp Penalty factors for each of the coefficients.
pen_function <- function(coefficients, alpha = 1.0, vp = 1.0) {
    sum(vp * (alpha * abs(coefficients) + (1-alpha)/2 * coefficients^2))
}

#' Elastic net deviance value
#'
#' Returns the elastic net deviance value.
#'
#' @param y Quantitative response variable.
#' @param mu Model's predictions for \code{y}.
#' @param weights Observation weights.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function.
dev_function <- function(y, mu, weights, family) {
    sum(family$dev.resids(y, mu, weights))
}


jerr.glmnetfit <- function (n, maxit, k = NULL) {
    if (n == 0) {
        list(n = 0, fatal = FALSE, msg = "")
    } else if (n > 0) {
        # fatal error
        fatal <- TRUE
        msg <- ifelse(n < 7777,
                      "Memory allocation error; contact package maintainer",
                      "Unknown error")
    } else {
        # non-fatal error
        fatal <- FALSE
        msg <- paste("Convergence for ", k, "th lambda value not reached after maxit=",
                     maxit, " iterations; solutions for larger lambdas returned",
                     sep = "")
    }
    list(n = n, fatal = fatal, msg = msg)
}
#' Get predictions from a \code{glmnetfit} fit object
#'
#' Gives fitted values, linear predictors, coefficients and number of non-zero
#' coefficients from a fitted \code{glmnetfit} object.
#'
#' @param object Fitted "glmnetfit" object.
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix. This argument is not used for \code{type =
#' c("coefficients","nonzero")}.
#' @param s Value(s) of the penalty parameter lambda at which predictions are
#' required. Default is the entire sequence used to create the model.
#' @param type Type of prediction required. Type "link" gives the linear
#' predictors (eta scale); Type "response" gives the fitted values (mu scale).
#' Type "coefficients" computes the coefficients at the requested values for s.
#' Type "nonzero" returns a list of the indices of the nonzero coefficients for
#' each value of s.
#' @param exact This argument is relevant only when predictions are made at values
#' of \code{s} (lambda) \emph{different} from those used in the fitting of the
#' original model. If \code{exact=FALSE} (default), then the predict function
#' uses linear interpolation to make predictions for values of \code{s} (lambda)
#' that do not coincide with those used in the fitting algorithm. While this is
#' often a good approximation, it can sometimes be a bit coarse. With
#' \code{exact=TRUE}, these different values of \code{s} are merged (and sorted)
#' with \code{object$lambda}, and the model is refit before predictions are made.
#' In this case, it is required to supply the original data x= and y= as additional
#' named arguments to predict() or coef(). The workhorse \code{predict.glmnet()}
#' needs to update the model, and so needs the data used to create it. The same
#' is true of weights, offset, penalty.factor, lower.limits, upper.limits if
#' these were used in the original call. Failure to do so will result in an error.
#' @param newoffset If an offset is used in the fit, then one must be supplied for
#' making predictions (except for type="coefficients" or type="nonzero").
#' @param ... This is the mechanism for passing arguments like \code{x=} when
#' \code{exact=TRUE}; see \code{exact} argument.
#'
#' @return The object returned depends on type.
#'
#' @method predict glmnetfit
#' @export
predict.glmnetfit <- function(object, newx, s = NULL,
                              type = c("link", "response", "coefficients", "nonzero"),
                              exact = FALSE, newoffset, ...) {
    type = match.arg(type)
    nfit <- NextMethod("predict")
    if (type == "response") {
        object$family$linkinv(nfit)
    } else {
        nfit
    }
}

#' Helper function to get etas (linear predictions)
#'
#' Given x, coefficients and intercept, return linear predictions. Wrapper that
#' works with both regular and sparse x. Only works for single set of coefficients
#' and intercept.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed to be standardized.
#' @param beta Feature coefficients.
#' @param a0 Intercept.
get_eta <- function(x, beta, a0) {
    if (inherits(x, "sparseMatrix")) {
        beta <- drop(beta)/attr(x, "xs")
        drop(x %*% beta - sum(beta * attr(x, "xm") ) + a0)
    } else {
        drop(x %*% beta + a0)
    }
}

#' Helper function to compute weighted mean and standard deviation
#'
#' Helper function to compute weighted mean and standard deviation.
#' Deals gracefully whether x is sparse matrix or not.
#'
#' @param x Observation matrix.
#' @param weights Optional weight vector.
#'
#' @return A list with components.
#' \item{mean}{vector of weighted means of columns of x}
#' \item{sd}{vector of weighted standard deviations of columns of x}
weighted_mean_sd <- function(x, weights=rep(1,nrow(x))){
    weights <- weights/sum(weights)
    xm <- drop(t(weights)%*%x)
    xv <- drop(t(weights)%*%scale(x,xm,FALSE)^2)
    xv[xv < 10*.Machine$double.eps] <- 0
    list(mean = xm, sd = sqrt(xv))
}

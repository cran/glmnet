#' Cox regression via penalized maximum likelihood using C++ engine
#'
#' This function calls the C++ implementation of Cox regression with
#' elastic net regularization. It handles both right-censored and
#' left-truncated (start, stop) survival data using the Breslow or Efron
#' method for ties. For stratified Cox models, it uses an IRLS approach
#' with integrated C++ gradient/Hessian computation.
#'
#' @param x Design matrix, of dimension nobs x nvars.
#' @param is.sparse Logical, is x a sparse matrix?
#' @param y Survival response variable, must be a Surv or stratifySurv object.
#' @param weights Observation weights.
#' @param offset Offset for the linear predictor.
#' @param alpha The elastic net mixing parameter.
#' @param nobs Number of observations.
#' @param nvars Number of variables.
#' @param jd Excluded variable indices (1-indexed, first element is count).
#' @param vp Penalty factors for each coefficient.
#' @param cl Coefficient limits matrix (2 x nvars).
#' @param ne Maximum number of variables in the model.
#' @param nx Maximum number of variables ever to be nonzero.
#' @param nlam Number of lambda values.
#' @param flmin Minimum lambda ratio.
#' @param ulam User-supplied lambda sequence.
#' @param thresh Convergence threshold.
#' @param isd Standardize flag.
#' @param vnames Variable names.
#' @param maxit Maximum number of iterations.
#' @param pb Progress bar object.
#' @param efron Logical; if TRUE use Efron method for ties, otherwise Breslow.
#'
#' @return An object of class "coxnet" with components:
#' \item{a0}{NULL (Cox model has no intercept)}
#' \item{beta}{Sparse coefficient matrix}
#' \item{df}{Number of nonzero coefficients per lambda}
#' \item{dim}{Dimension of coefficient matrix}
#' \item{lambda}{Lambda sequence used}
#' \item{dev.ratio}{Fraction of null deviance explained}
#' \item{nulldev}{Null deviance}
#' \item{npasses}{Number of coordinate descent passes}
#' \item{jerr}{Error code}
#' \item{offset}{Logical indicating if offset was used}
#'
#' @keywords internal
coxnet <- function(x, is.sparse, y, weights, offset, alpha, nobs, nvars, jd, vp,
                   cl, ne, nx, nlam, flmin, ulam, thresh, isd, vnames, maxit, pb,
                   efron = FALSE) {

    # Preserve strata attribute before parsing
    strata <- NULL
    if ("strata" %in% names(attributes(y))) {
        strata <- attr(y, "strata")
    }

    # Parse survival response
    y <- response.coxnet(y)
    if (!is.matrix(y))
        stop("Cox model requires a matrix response; a 'Surv' object suffices",
             call. = FALSE)

    # Create strata vector for C++ (empty for non-stratified to save memory)
    if (is.null(strata)) {
        strata <- integer(0)
    } else {
        # Convert to integer factor levels (1-indexed)
        strata <- as.integer(as.factor(strata))
    }

    # Unified C++ path for all Cox models (stratified and non-stratified)
    # Determine if this is (start, stop) or right-censored data
    if (ncol(y) == 2) {
        # Right-censored: y has columns "time" and "status"
        if (!all(match(c("time", "status"), dimnames(y)[[2]], 0)))
            stop("Cox model requires columns 'time' (>0) and 'status' (binary)",
                 call. = FALSE)
        start_time <- rep(0.0, nobs)
        stop_time <- as.double(y[, "time"])
        status <- as.integer(y[, "status"])
    } else if (ncol(y) == 3) {
        # Left-truncated (start, stop]: y has columns "start", "stop", "status"
        if (!all(match(c("start", "stop", "status"), dimnames(y)[[2]], 0)))
            stop("Left-truncated Cox model requires columns 'start', 'stop', and 'status'",
                 call. = FALSE)
        start_time <- as.double(y[, "start"])
        stop_time <- as.double(y[, "stop"])
        status <- as.integer(y[, "status"])
    } else {
        stop("Unexpected survival response format", call. = FALSE)
    }

    # Validate times
    if (any(stop_time <= 0))
        stop("Non-positive event times encountered; not permitted for Cox family",
             call. = FALSE)
    if (any(start_time < 0))
        stop("Negative start times encountered; not permitted for Cox family",
             call. = FALSE)
    if (any(start_time >= stop_time))
        stop("Start time must be less than stop time", call. = FALSE)

    # Ensure proper types
    maxit <- as.integer(maxit)
    weights <- as.double(weights)

    # Handle offset
    if (is.null(offset)) {
        offset <- rep(0.0, nobs)
        is.offset <- FALSE
    } else {
        storage.mode(offset) <- "double"
        is.offset <- TRUE
    }

    # Call unified C++ implementation (handles both stratified and non-stratified)
    if (is.sparse) {
        fit <- spcoxnet_exp(
            parm = alpha,
            x = x,
            start = start_time,
            stop = stop_time,
            status = status,
            strata = strata,
            efron = efron,
            g = offset,
            w = weights,
            jd = jd,
            vp = vp,
            cl = cl,
            ne = ne,
            nx = nx,
            nlam = nlam,
            flmin = flmin,
            ulam = ulam,
            thr = thresh,
            isd = isd,
            maxit = maxit,
            pb = pb,
            lmu = integer(1),
            a0 = double(nlam),
            ca = matrix(0.0, nx, nlam),
            ia = integer(nx),
            nin = integer(nlam),
            nulldev = double(1),
            dev = double(nlam),
            alm = double(nlam),
            nlp = integer(1),
            jerr = integer(1)
        )
    } else {
        fit <- coxnet_exp(
            parm = alpha,
            x = x,
            start = start_time,
            stop = stop_time,
            status = status,
            strata = strata,
            efron = efron,
            g = offset,
            w = weights,
            jd = jd,
            vp = vp,
            cl = cl,
            ne = ne,
            nx = nx,
            nlam = nlam,
            flmin = flmin,
            ulam = ulam,
            thr = thresh,
            isd = isd,
            maxit = maxit,
            pb = pb,
            lmu = integer(1),
            a0 = double(nlam),
            ca = matrix(0.0, nx, nlam),
            ia = integer(nx),
            nin = integer(nlam),
            nulldev = double(1),
            dev = double(nlam),
            alm = double(nlam),
            nlp = integer(1),
            jerr = integer(1)
        )
    }

    # Handle errors
    if (fit$jerr != 0) {
        errmsg <- jerr(fit$jerr, maxit, pmax = nx, family = "cox")
        if (errmsg$fatal)
            stop(errmsg$msg, call. = FALSE)
        else
            warning(errmsg$msg, call. = FALSE)
    }

    # Format output
    # Cox model has no intercept, so set a0 to NULL before calling getcoef
    fit$a0 <- NULL
    outlist <- getcoef(fit, nvars, nx, vnames)
    dev <- fit$dev[seq(fit$lmu)]
    outlist <- c(outlist, list(
        dev.ratio = dev,
        nulldev = fit$nulldev,
        npasses = fit$nlp,
        jerr = fit$jerr,
        offset = is.offset
    ))
    class(outlist) <- "coxnet"
    outlist
}

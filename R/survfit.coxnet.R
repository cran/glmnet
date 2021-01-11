#' Compute a survival curve from a coxnet object
#'
#' Computes the predicted survivor function for a Cox proportional hazards
#' model with elastic net penalty.
#'
#' To be consistent with other functions in \code{glmnet}, if \code{s}
#' is not specified, survival curves are returned for the entire lambda
#' sequence. This is not recommended usage: it is best to call
#' \code{survfit.coxnet} with a single value of the penalty parameter
#' for the \code{s} option.
#'
#' @param formula A class \code{coxnet} object.
#' @param s Value(s) of the penalty parameter lambda at which the survival
#' curve is required. Default is the entire sequence used to create the model.
#' However, it is recommended that \code{survfit.coxnet} is called for
#' a single penalty parameter.
#' @param ... This is the mechanism for passing additional arguments like
#' (i) x= and y= for the x and y used to fit the model,
#' (ii) weights= and offset= when the model was fit with these options,
#' (iii) arguments for new data (newx, newoffset, newstrata), and
#' (iv) arguments to be passed to survfit.coxph().
#'
#' @return If \code{s} is a single value, an object of class "survfitcox"
#' and "survfit" containing one or more survival curves. Otherwise, a list
#' of such objects, one element for each value in \code{s}.
#' Methods defined for survfit objects are print, summary and plot.
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
#' y <- survival::Surv(ty, tcens)
#' fit1 <- glmnet(x, y, family = "cox")
#'
#' # survfit object for Cox model where lambda = 0.1
#' sf1 <- survival::survfit(fit1, s = 0.1, x = x, y = y)
#' plot(sf1)
#'
#' # example with new data
#' sf2 <- survival::survfit(fit1, s = 0.1, x = x, y = y, newx = x[1:3, ])
#' plot(sf2)
#'
#' # example with strata
#' y2 <- stratifySurv(y, rep(1:2, length.out = nobs))
#' fit2 <- glmnet(x, y2, family = "cox")
#' sf3 <- survival::survfit(fit2, s = 0.1, x = x, y = y2)
#' sf4 <- survival::survfit(fit2, s = 0.1, x = x, y = y2,
#'                newx = x[1:3, ], newstrata = c(1, 1, 1))
#'
#' @importFrom survival coxph survfit
#' @method survfit coxnet
#' @export
survfit.coxnet <- function(formula, s = NULL, ...) {
  this.call <- match.call()
  object <- formula
  args <- list(...)

  if (!("x" %in% names(args)) || !("y" %in% names(args)))
    stop(paste0("the design matrix x and response y used to fit the model ",
                "need to be passed"))
  y <- args$y

  # if s is NULL, get the whole lambda sequence
  if (is.null(s)) s <- object$lambda

  # check that the required arguments are passed for running coxph
  check_dots(object, ...,
             need = c("offset", "weights"),
             error_start = "used survfit.coxnet()",
             error_end = "in order to run survfit.coxnet")
  # if user wants to run survfit on new data, make sure the arguments required
  # for making new predictions are provided
  if ("newx" %in% names(args)) {
    check_dots(object, ...,
               need = c("offset"),
               error_start = "used survfit.coxnet() with newx argument",
               prefix = "new",
               error_end = "in order to predict on new data")

    if ("strata" %in% names(attributes(y)) &&
        !("newstrata" %in% names(args)))
      stop(paste0("used survfit.coxnet() with newx argument and ",
                  "stratified Cox model was fit, need newstrata argument ",
                  "in order to predict on new data"))
  }

  res <- list()
  for (i in seq_along(s)) {
    # "hack": run coxph with 0 iterations
    coxphmod <- mycoxph(object, s = s[i], ...)

    # If newx is provided, we need to compute the predections at these
    # observations and pass it as an argument to the eventual survfit call
    # as `newdata`. In addition, we have to handle additional options that
    # could have been passed: offset
    if ("newx" %in% names(args)) {
      current_args <- mycoxpred(object, s = s[i], ...)
    } else {
      current_args <- args
    }

    # make the call to survfit.coxph
    current_args$formula <- coxphmod
    current_args$se.fit <- FALSE  # doesn't make sense to compute SEs
    sf <- do.call(survfit, current_args)
    sf$call <- this.call
    res[[i]] <- sf
  }
  if (length(s) > 1) {
    return(res)
  } else {
    return(res[[1]])
  }
}

#' Helper function to fit coxph model for survfit.coxnet
#'
#' This function constructs the coxph call needed to run the "hack" of
#' coxph with 0 iterations. It's a separate function as we have to deal with
#' function options like strata, offset and observation weights.
#'
#' @param object A class \code{coxnet} object.
#' @param s The value of the penalty parameter lambda at which the survival
#' curve is required.
#' @param ... The same ... that was passed to survfit.coxnet.
#'
#' @importFrom survival coxph strata
#' @importFrom stats as.formula
mycoxph <- function(object, s, ...) {
  args <- list(...)
  x <- args$x
  y <- args$y
  glmnet_call_names <- names(object$call)[-1]

  # predict from coxnet model for original data frame x with s, gives fitted
  # linear predictor
  args$object <- object
  args$newx <- x
  args$s <- s
  if ("offset" %in% glmnet_call_names) {
    args$newoffset <- rep(0, length.out = nrow(x))
  }
  eta <- do.call(predict, args)

  # construct list of arguments for coxph() call based on which special
  # arguments were used in the original glmnet() call
  coxphargs <- list(formula = "y ~ X1", data = data.frame(y, eta),
                    init = 1, iter = 0)

  if ("strata" %in% names(attributes(y))) {
    coxphargs$data$strata <- attr(y, "strata")
    coxphargs$formula <- paste(coxphargs$formula, "+ strata(strata)")
  }
  if ("weights" %in% glmnet_call_names) {
    coxphargs$weights <- args$weights
  }
  if ("offset" %in% glmnet_call_names) {
    coxphargs$data$offset <- args$offset
    coxphargs$formula <- paste(coxphargs$formula, "+ offset(offset)")
  }
  coxphargs$formula <- as.formula(coxphargs$formula)
  coxphmod <- do.call(coxph, coxphargs)
  return(coxphmod)
}

#' Helper function to amend ... for new data in survfit.coxnet
#'
#' This function amends the function arguments passed to survfit.coxnet
#' via ... if new data was passed to survfit.coxnet. It's a separate
#' function as we have to deal with function options like newstrata
#' and newoffset.
#'
#' @param object A class \code{coxnet} object.
#' @param s The response for the fitted model.
#' @param ... The same ... that was passed to survfit.coxnet.
#'
#' @importFrom stats predict
mycoxpred <- function(object, s, ...) {
  args <- list(...)

  if ("newoffset" %in% names(args)) {
    new_eta <- predict(object, newx = args$newx, s = s,
                       newoffset = rep(0, nrow(args$newx)))
  } else {
    new_eta <- predict(object, newx = args$newx, s = s)
  }
  new_df <- data.frame(new_eta)
  if ("newoffset" %in% names(args)) new_df$offset <- args$newoffset

  if ("newstrata" %in% names(args)) {
    new_df$strata <- args$newstrata
    args$newstrata <- NULL
  }

  args$newdata <- new_df
  return(args)
}

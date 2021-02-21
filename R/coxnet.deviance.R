#' Compute deviance for Cox model
#' 
#' Compute the deviance (-2 log partial likelihood) for Cox model.
#'
#' Computes the deviance for a single set of predictions, or for a matrix 
#' of predictions. The user can either supply the predictions
#' directly through the \code{pred} option, or by supplying the \code{x} matrix
#' and \code{beta} coefficients. Uses the Breslow approach to ties.
#' 
#' The function first checks if \code{pred} is passed: if so, it is used as
#' the predictions. If \code{pred} is not passed but \code{x} and \code{beta}
#' are passed, then these values are used to compute the predictions. If 
#' neither \code{x} nor \code{beta} are passed, then the predictions are all
#' taken to be 0.
#' 
#' \code{coxnet.deviance()} is a wrapper: it calls the appropriate internal
#' routine based on whether the response is right-censored data or 
#' (start, stop] survival data. 
#' 
#' @aliases coxnet.deviance
#' @param pred Fit vector or matrix (usually from glmnet at a particular 
#' lambda or a sequence of lambdas).
#' @param y Survival response variable, must be a \code{Surv} or 
#' \code{stratifySurv} object.
#' @param x Optional \code{x} matrix, to be supplied if \code{pred = NULL}.
#' @param offset Optional offset vector.
#' @param weights Observation weights (default is all equal to 1).
#' @param std.weights If TRUE (default), observation weights are standardized
#' to sum to 1.
#' @param beta Optional coefficient vector/matrix, to be supplied if 
#' \code{pred = NULL}.
#' 
#' @return A vector of deviances, one for each column of predictions.
#' 
#' @examples
#' set.seed(1)
#' eta <- rnorm(10)
#' time <- runif(10, min = 1, max = 10)
#' d <- ifelse(rnorm(10) > 0, 1, 0)
#' y <- survival::Surv(time, d)
#' coxnet.deviance(pred = eta, y = y)
#' 
#' # if pred not provided, it is set to zero vector
#' coxnet.deviance(y = y)
#' 
#' # example with x and beta
#' x <- matrix(rnorm(10 * 3), nrow = 10)
#' beta <- matrix(1:3, ncol = 1)
#' coxnet.deviance(y = y, x = x, beta = beta)
#' 
#' # example with (start, stop] data
#' y2 <- survival::Surv(time, time + runif(10), d)
#' coxnet.deviance(pred = eta, y = y2)
#' 
#' # example with strata
#' y2 <- stratifySurv(y, rep(1:2, length.out = 10))
#' coxnet.deviance(pred = eta, y = y2)
#' 
#' @seealso \code{coxgrad}
#' @keywords Cox model
#'
#' @export
coxnet.deviance <- function(pred = NULL, y, x = NULL, offset = NULL, 
                            weights = NULL, std.weights = TRUE, beta = NULL) {
  y <- response.coxnet(y)
  
  # if y has 2 columns, it is right-censored data
  # if y has 3 columns, it is (start, stop] data
  # otherwise, throw error
  if (ncol(y) == 2) {
    return(coxnet.deviance0(pred = pred, y = y, x = x, offset = offset,
                            weights = weights, std.weights = std.weights,
                            beta = beta))
  } else if (ncol(y) == 3) {
    return(coxnet.deviance3(pred = pred, y = y, x = x, offset = offset, 
                            weights = weights, std.weights = std.weights,
                            beta = beta))
  } else {
    stop("Response y should have 2 or 3 columns")
  }
}

# coxnet.deviance routine for right-censored data
coxnet.deviance0 <- function(pred = NULL, y, x = NULL, offset = NULL, 
                             weights = NULL, std.weights = TRUE, beta = NULL) {
  ty <- y[, "time"]
  tevent <- y[, "status"]
  ty <- ty + (1 - tevent) * 100 * .Machine$double.eps
  nobs <- as.integer(length(ty))
  
  # hack for the case where user passes in x as sparse matrix
  if (!is.null(x) && inherits(x, "sparseMatrix")) {
    if (is.null(beta))
      stop("if x is passed, beta must also be passed")
    pred <- as.matrix(x %*% beta)
    return(coxnet.deviance0(pred = pred, y = y, offset = offset,
                            weights = weights, std.weights = std.weights))
  }
  
  # Sort out the pred, x and beta options.
  # If user provided `pred`, we let x = pred and beta = identity matrix.
  # This allows us to use the loglike Fortran routine to compute the
  # partial log likelihood.
  # In the end, only x and beta are passed to the Fortran routine.
  if (!is.null(pred)) {
    x <- as.matrix(pred)
    nvec <- ncol(x)
    beta <- diag(nvec)
    nvars <- as.integer(nvec)
  } else if (is.null(x) && is.null(beta)) {
    x <- matrix(0, nrow = nobs, ncol = 1)
    beta <- double(0)
    nvec <- 1
    nvars <- as.integer(0)
  } else if (!is.null(x) && !is.null(beta)) {
    x <- as.matrix(x)
    beta <- as.matrix(beta)
    nvec <- ncol(beta)
    nvars <- nrow(beta)
  } else {
    stop("user must pass either `pred`, or both `x` and `beta`")
  }
  storage.mode(x) <- "double"
  storage.mode(beta) <- "double"
  nvec <- as.integer(nvec)
  nvars <- as.integer(nvars)
  
  # normalize weights to sum to nobs
  if (is.null(weights))
    weights <- rep(1, nobs)
  else {
    if (std.weights) weights <- nobs * weights / sum(weights)
    weights <- as.double(weights)
  }
  
  if (is.null(offset))
    offset <- rep(0, nobs)
  else offset <- as.double(offset)
  
  # extract strata (if any)
  if ("strata" %in% names(attributes(y))) {
    strata <- attr(y, "strata")
  } else {
    strata <- rep(1, nobs)
  }
  if (length(strata) != nobs) stop("length of strata != nobs")
  
  # if all in same strata, do the deviance computation
  # if not, take the sum of the strata-level deviances
  if (length(unique(strata)) == 1) {
    ### Compute saturated loglikelihood
    wd <- weights[tevent == 1]
    tyd <- ty[tevent == 1]
    if (any(duplicated(tyd))) {
      wd <- tapply(wd, tyd, sum)
    }
    wd <- wd[wd > 0]
    lsat <- -sum(wd * log(wd))
    ####
    
    fit <- .Fortran("loglike", nobs, nvars, x, ty, tevent, offset,
                    weights, nvec, beta, flog = double(nvec), jerr = integer(1),
                    PACKAGE = "glmnet")
    if (fit$jerr != 0) {
      errmsg <- jerr(fit$jerr, maxit = 0, pmax = 0, family = "cox")
      if (errmsg$fatal)
        stop(errmsg$msg, call. = FALSE)
      else warning(errmsg$msg, call. = FALSE)
    }
    return(2 * (lsat - fit$flog))
  } else {
    # more than one strata provided: return the sum of strata-level deviances
    tot_dev <- 0
    for (i in unique(strata)) {
      ii <- which(strata == i)
      tot_dev <- tot_dev + 
        coxnet.deviance0(y = y[ii, , drop = FALSE], x = x[ii, , drop = FALSE], 
                         beta = beta, offset = offset[ii], 
                         weights = weights[ii], std.weights = FALSE)
    }
    return(tot_dev)
  }
}

# coxnet.deviance2 gives the same output as coxnet.deviance0()
# but is written completely in R. It is not called by 
# coxnet.deviance(), and is kept in the package for completeness.
coxnet.deviance2 <- function(pred = NULL, y, x = NULL, offset = NULL, 
                             weights = NULL, std.weights = TRUE, beta = NULL) {
  if (!is.Surv(y)) stop("y must be a Surv object")
  nobs <- nrow(y)
  
  # if pred is NULL, use beta and x to compute pred
  # if beta is NULL too, set pred to all zeros
  if (is.null(pred)) {
    if ((!is.null(x) && is.null(beta)) || (is.null(x) && !is.null(beta)))
      stop("user must pass either `pred`, or both `x` and `beta`")
    if (is.null(beta)) {
      pred <- rep(0, times = nobs)
    } else {
      pred <- x %*% beta
    }
  }
  
  # if more than one column of predictions is passed, run coxnet.deviance2()
  # for each column
  if (!is.null(ncol(pred)) && ncol(pred) > 1) {
    return(sapply(seq(ncol(pred)),
                  function(j) coxnet.deviance2(
                    pred = pred[, j], y = y, offset = offset, 
                    weights = weights, std.weights = std.weights)))
  } else {
    # check that pred is of the right length
    if(length(pred) != nobs) stop("pred and y must have the same length")
    
    # normalize weights to sum to nobs
    if (is.null(weights))
      w <- rep(1, nobs)
    else {
      if (length(weights) != nobs) stop("weights and y must have the same length")
      if (std.weights) {
        w <- nobs * weights / sum(weights)
      } else {
        w <- weights
      }
    }
    
    # if there's an offset, add it to the pred vector
    if (is.null(offset)) {
      offset <- rep(0, nobs)
    } else {
      if (length(offset) != nobs) stop("offset and y must have the same length")
      pred <- pred + offset
    }
    
    # extract strata (if any)
    if ("strata" %in% names(attributes(y))) {
      strata <- attr(y, "strata")
    } else {
      strata <- rep(1, nobs)
    }
    if (length(strata) != nobs) stop("length of strata != nobs")
    
    # if all in same strata, do the deviance computation
    # if not, take the sum of the strata-level deviances
    if (length(unique(strata)) == 1) {
      time <- y[, "time"]
      d <- y[, "status"]
      
      ### Compute saturated loglikelihood
      wd <- w[d == 1]
      tyd <- time[d == 1]
      if (any(duplicated(tyd))) {
        wd <- tapply(wd, tyd, sum)
      }
      wd <- wd[wd > 0]
      lsat <- -sum(wd * log(wd))
      ####
      
      # order time, d, pred and w in ascending time order
      # for tied times, all deaths come before censored observations
      if ("stop_time" %in% names(attributes(y))) {
        o <- attr(y, "stop_time")
      } else {
        o <- order(time, d, decreasing = c(FALSE, TRUE))
      }
      time <- time[o]
      d <- d[o]
      pred <- pred[o]
      w <- w[o]
      
      ### See if there are dups in death times
      dups <- fid(time[d==1],seq(length(d))[d==1])
      dd <- d
      ww <- w
      
      ### next code replaces each sequence of tied death indicators by a new
      ### sequence where only the first is a 1 and the rest are zero. This
      ### makes the accounting in the following step work properly we also
      ### sums the weights in each of the tied death sets, and assign that
      ### weight to the first
      if(!is.null(ties<-dups$index_ties)){
        dd[unlist(ties)]=0
        dd[dups$index_first]=1
        wsum=sapply(ties,function(i,w)sum(w[i]),ww)
        tie1=sapply(ties,function(i)i[1])
        ww[tie1]=wsum
      }
      
      # compute the sum inside the log term of the partial likelihood
      w_exp_pred <- w * exp(pred)
      rsk <- rev(cumsum(rev(w_exp_pred)))
      
      # take just the terms related to actual death times
      log_terms <- (ww * log(rsk))[dd > 0]
      loglik <- sum((w * pred)[d > 0]) - sum(log_terms)
      
      return(2 * (lsat -loglik))
    } else {
      # more than one strata provided: return the sum of strata-level 
      # deviances
      tot_dev <- 0
      for (i in unique(strata)) {
        ii <- which(strata == i)
        tot_dev <- tot_dev + 
          coxnet.deviance2(pred = pred[ii], y = y[ii, , drop = FALSE], 
                           offset = NULL, weights = w[ii], std.weights = FALSE)
      }
      return(tot_dev)
    }
  }
}

# coxnet.deviance routine for (start, stop] data
coxnet.deviance3 <- function(pred = NULL, y, x = NULL, offset = NULL, 
                             weights = NULL, std.weights = TRUE, beta = NULL) {
  if (!is.Surv(y)) stop("y must be a Surv object")
  nobs <- nrow(y)
  
  # if pred is NULL, use beta and x to compute pred
  # if beta is NULL too, set pred to all zeros
  if (is.null(pred)) {
    if ((!is.null(x) && is.null(beta)) || (is.null(x) && !is.null(beta)))
      stop("user must pass either `pred`, or both `x` and `beta`")
    if (is.null(beta)) {
      pred <- rep(0, times = nobs)
    } else {
      pred <- x %*% beta
    }
  }
  
  # if more than one column of predictions is passed, run coxnet.deviance3()
  # for each column
  if (!is.null(ncol(pred)) && ncol(pred) > 1) {
    return(sapply(seq(ncol(pred)),
                  function(j) coxnet.deviance3(
                    pred = pred[, j], y = y, offset = offset, 
                    weights = weights, std.weights = std.weights)))
  } else {
    # check that pred is of the right length
    if(length(pred) != nobs) stop("pred and y must have the same length")
    
    # normalize weights to sum to nobs
    if (is.null(weights))
      w <- rep(1, nobs)
    else {
      if (length(weights) != nobs) stop("weights and y must have the same length")
      if (std.weights) {
        w <- nobs * weights / sum(weights)
      } else {
        w <- weights
      }
    }
    
    # if there's an offset, add it to the pred vector
    if (is.null(offset)) {
      offset <- rep(0, nobs)
    } else {
      if (length(offset) != nobs) stop("offset and y must have the same length")
      pred <- pred + offset
    }
    
    # extract strata (if any)
    if ("strata" %in% names(attributes(y))) {
      strata <- attr(y, "strata")
    } else {
      strata <- rep(1, nobs)
    }
    if (length(strata) != nobs) stop("length of strata != nobs")
    
    # if all in same strata, do the deviance computation
    # if not, take the sum of the strata-level deviances
    if (length(unique(strata)) == 1) {
      start_time <- y[, "start"]
      stop_time <- y[, "stop"]
      d <- y[, "status"]
      
      ### Compute saturated loglikelihood
      wd <- w[d == 1]
      tyd <- stop_time[d == 1]
      if (any(duplicated(tyd))) {
        wd <- tapply(wd, tyd, sum)
      }
      wd <- wd[wd > 0]
      lsat <- -sum(wd * log(wd))
      ####
      
      # get ordering for stop time (ascending, deaths before censored), and
      # start time (ascending)
      if ("stop_time" %in% names(attributes(y))) {
        stop_o <- attr(y, "stop_time")
      } else {
        stop_o <- order(stop_time, d, decreasing = c(FALSE, TRUE))
      }
      if ("start_time" %in% names(attributes(y))) {
        start_o <- attr(y, "start_time")
      } else {
        start_o <- order(start_time, decreasing = c(FALSE))
      }
      
      # keep a set of values which are ordered by start time
      w_exp_pred_start <- (w * exp(pred))[start_o]
      start_time_start <- start_time[start_o]
      
      # reorder everything by stop time
      start_time <- start_time[stop_o]
      stop_time  <- stop_time[stop_o]
      d          <- d[stop_o]
      pred       <- pred[stop_o]
      w          <- w[stop_o]
      
      ### See if there are dups in death times
      dups <- fid(stop_time[d == 1], seq(length(d))[d == 1])
      dd <- d
      ww <- w
      
      ### next code replaces each sequence of tied death indicators by a new
      ### sequence where only the first is a 1 and the rest are zero. This
      ### makes the accounting in the following step work properly we also
      ### sums the weights in each of the tied death sets, and assign that
      ### weight to the first
      if(!is.null(ties<-dups$index_ties)){
        dd[unlist(ties)]=0
        dd[dups$index_first]=1
        wsum=sapply(ties,function(i,w)sum(w[i]),ww)
        tie1=sapply(ties,function(i)i[1])
        ww[tie1]=wsum
      }
      
      # compute risk set sums rsk[i] = \sum_{j in R_i} w_j exp(eta_j)
      # where i indexes the observations. (In the end, we will only care
      # about the indices i which have actual death times.)
      rsk <- rev(cumsum(rev(w * exp(pred))))
      current_sum <- 0
      stop_idx <- nobs; start_idx <- nobs
      while (stop_idx > 0 && start_idx > 0) {
        if (start_time_start[start_idx] < stop_time[stop_idx]) {
          # current start time belongs in risk set ending in stop time,
          # so we should remove the current cumulative sum and consider
          # the next risk set
          rsk[stop_idx] <- rsk[stop_idx] - current_sum
          stop_idx <- stop_idx - 1
        } else {
          # current start time does not belong in risk set ending in stop
          # time, so we should add it to current_sum and check if the
          # start time before it should also be added
          current_sum <- current_sum + w_exp_pred_start[start_idx]
          start_idx <- start_idx - 1
        }
      }
      
      log_terms <- ww[dups$index_first] * (log(rsk[dd == 1]))
      loglik <- sum((w * pred)[d > 0]) - sum(log_terms)
      return(2 * (lsat -loglik))
    } else {
      # more than one strata provided: return the sum of strata-level 
      # deviances
      tot_dev <- 0
      for (i in unique(strata)) {
        ii <- which(strata == i)
        tot_dev <- tot_dev + 
          coxnet.deviance3(pred = pred[ii], y = y[ii, , drop = FALSE], 
                           offset = NULL, weights = w[ii], std.weights = FALSE)
      }
      return(tot_dev)
    }
  }
}
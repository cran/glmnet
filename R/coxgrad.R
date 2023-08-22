#' Compute gradient for Cox model
#'
#' Compute the gradient of the log partial likelihood at a particular fit for Cox
#' model.
#'
#' Compute a gradient vector at the fitted vector for the log partial likelihood.
#' This is like a residual vector, and useful for manual screening of
#' predictors for \code{glmnet} in applications where \code{p} is very large
#' (as in GWAS). Uses the Breslow approach to ties.
#'
#' This function is essentially a wrapper: it checks whether the response
#' provided is right-censored or (start, stop] survival data, and calls the
#' appropriate internal routine.
#'
#' @aliases coxgrad
#' @param eta Fit vector (usually from glmnet at a particular lambda).
#' @param y Survival response variable, must be a \code{Surv} or
#' \code{stratifySurv} object.
#' @param w Observation weights (default is all equal to 1).
#' @param std.weights If TRUE (default), observation weights are standardized
#' to sum to 1.
#' @param diag.hessian If \code{TRUE}, compute the diagonal of the Hessian
#' of the log partial likelihood as well. Default is \code{FALSE}.
#'
#' @return A single gradient vector the same length as \code{eta}. If
#' \code{diag.hessian=TRUE}, the diagonal of the Hessian is
#' included as an attribute "diag_hessian".
#'
#' @examples
#' set.seed(1)
#' eta <- rnorm(10)
#' time <- runif(10, min = 1, max = 10)
#' d <- ifelse(rnorm(10) > 0, 1, 0)
#' y <- survival::Surv(time, d)
#' coxgrad(eta, y)
#'
#' # return diagonal of Hessian as well
#' coxgrad(eta, y, diag.hessian = TRUE)
#'
#' # example with (start, stop] data
#' y2 <- survival::Surv(time, time + runif(10), d)
#' coxgrad(eta, y2)
#'
#' # example with strata
#' y2 <- stratifySurv(y, rep(1:2, length.out = 10))
#' coxgrad(eta, y2)
#'
#' @seealso \code{coxnet.deviance}
#' @keywords Cox model
#'
#' @export
coxgrad <- function(eta, y, w, std.weights = TRUE, diag.hessian = FALSE) {
    # if y has 2 columns, it is right-censored data
    # if y has 3 columns, it is (start, stop] data
    # otherwise, throw errors
    if (ncol(y) == 2) {
        return(coxgrad2(eta, y, w, std.weights, diag.hessian))
    } else if (ncol(y) == 3) {
        return(coxgrad3(eta, y, w, std.weights, diag.hessian))
    } else {
        stop("Response y should have 2 or 3 columns")
    }
}

# coxgrad routine for right-censored data
coxgrad2 <- function(eta, y, w, std.weights = TRUE, diag.hessian = FALSE) {
    if (missing(w)) w=rep(1,length(eta))
    if (std.weights) w=w/sum(w)
    nobs <- nrow(y)

    # extract strata (if any)
    if ("strata" %in% names(attributes(y))) {
        strata <- attr(y, "strata")
    } else {
        strata <- rep(1, nobs)
    }
    if (length(strata) != nobs) stop("length of strata != nobs")

    # if all in same strata, do the computations
    # if not, do strata-level computations and concatenate
    if (length(unique(strata)) == 1) {
        time <- y[, "time"]
        d    <- y[, "status"]
        eta <- scale(eta, TRUE, FALSE)  # center eta so exponents are not too large

        # order exp(eta), time, d and w in ascending time order
        # for tied times, all deaths come before censored observations
        if ("stop_time" %in% names(attributes(y))) {
            o <- attr(y, "stop_time")
        } else {
            o <- order(time, d, decreasing = c(FALSE, TRUE))
        }
        exp_eta <- exp(eta)[o]
        time <- time[o]
        d <- d[o]
        w <- w[o]
        rskden <- rev(cumsum(rev(exp_eta*w))) ##reverse order inside;last guy is in all the risk sets

        ### See if there are dups in death times
        dups <- fid(time[d == 1],seq(length(d))[d == 1])
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

        ### Get counts over risk sets at each death time
        rskcount=cumsum(dd)#this says how many of the risk sets each observation is in; 0 is none
        ### We now form partial sums of the 1/den just at the risk sets
        rskdeninv=cumsum((ww/rskden)[dd==1])
        ### pad with a zero, so we can index it
        rskdeninv=c(0,rskdeninv)

        ### compute gradient for each obs
        grad <- w * (d - exp_eta * rskdeninv[rskcount+1])
        grad[o] <- grad

        # if diag.hessian = TRUE, return the diagonal of the hessian too
        if (diag.hessian) {
            rskdeninv2 <- cumsum((ww/(rskden^2))[dd==1])
            rskdeninv2 <- c(0, rskdeninv2)
            w_exp_eta <- w * exp_eta
            diag_hessian <- w_exp_eta^2 * rskdeninv2[rskcount+1] - w_exp_eta * rskdeninv[rskcount+1]
            diag_hessian[o] <- diag_hessian
            attr(grad, "diag_hessian") <- diag_hessian
        }
        return(grad)
    } else {
        # more than one strata provided: compute strata-level values and
        # concatenate
        overall_grad <- rep(NA, nobs)
        if (diag.hessian) overall_diag_hessian <- rep(NA, nobs)
        for (i in unique(strata)) {
            ii <- which(strata == i)
            strata_res <- coxgrad2(eta[ii], y[ii, , drop = FALSE], w[ii],
                                   std.weights = FALSE, diag.hessian = diag.hessian)
            overall_grad[ii] <- strata_res
            if (diag.hessian) {
                overall_diag_hessian[ii] <- attr(strata_res, "diag_hessian")
            }
        }
        if (diag.hessian) {
            attr(overall_grad, "diag_hessian") <- overall_diag_hessian
        }
        return(overall_grad)
    }
}

# coxgrad routine for (start, stop] data
coxgrad3 <- function(eta, y, w, std.weights = TRUE, diag.hessian = FALSE) {
    if (missing(w)) w=rep(1,length(eta))
    if (std.weights) w=w/sum(w)
    nobs <- nrow(y)

    # extract strata (if any)
    if ("strata" %in% names(attributes(y))) {
        strata <- attr(y, "strata")
    } else {
        strata <- rep(1, nobs)
    }
    if (length(strata) != nobs) stop("length of strata != nobs")

    # if all in same strata, do the computations
    # if not, do strata-level computations and concatenate
    if (length(unique(strata)) == 1) {
        start_time <- y[, "start"]
        stop_time  <- y[, "stop"]
        d          <- y[, "status"]
        eta <- scale(eta, TRUE, FALSE)  # center eta so exponents are not too large

        # get ordering for stop time (ascending, deaths before censored),
        # start time (ascending), and match info if cached
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
        if ("ss_match" %in% names(attributes(y))) {
            ss_match <- attr(y, "ss_match")
        } else {
            ss_match <- match(start_o, stop_o)
        }
    ## set weights to zero for any observation with start time greater than largest death time
    last_death <- max(stop_time[d==1])
    w[start_time >= last_death] <- 0
    ##

        # keep a set of values which are ordered by start time
        w_exp_eta_start <- (w * exp(eta))[start_o]
        start_time_start <- start_time[start_o]

        # reorder everything by stop time
        exp_eta <- exp(eta)[stop_o]
        start_time <- start_time[stop_o]
        stop_time <- stop_time[stop_o]
        d <- d[stop_o]
        w <- w[stop_o]

        ### See if there are dups in death times
        dups <- fid(stop_time[d == 1],seq(length(d))[d == 1])
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

        # compute risk set sums rskden[i] = \sum_{j in R_i} w_j exp(eta_j)
        # where i indexes the observations. (In the end, we will only care
        # about the indices i which have actual death times.)
        rskden <- rev(cumsum(rev(exp_eta*w)))
        current_sum <- 0
        death_time <- stop_time[dups$index_first]
        ndeaths <- length(death_time)
        death_idx <- ndeaths; start_idx <- nobs
        while (death_idx > 0 && start_idx > 0) {
            if (start_time_start[start_idx] < death_time[death_idx]) {
                # current start time belongs in risk set ending in stop time,
                # so we should remove the current cumulative sum and consider
                # the next risk set
                stop_idx <- dups$index_first[death_idx]
                rskden[stop_idx] <- rskden[stop_idx] - current_sum
                death_idx <- death_idx - 1
            } else {
                # current start time does not belong in risk set ending in stop
                # time, so we should add it to current_sum and check if the
                # start time before it should also be added
                current_sum <- current_sum + w_exp_eta_start[start_idx]
                start_idx <- start_idx - 1
            }
        }

        # compute the terms rskterm[k] = \sum_{i in C_k} d[i] / rskden[i] and
        # rskterm2[k] = \sum_{i in C_k} d[i] / rskden[i]^2.
        # Here, k indexes the observations, index i runs over the unique death
        # times.
        rskfactor  <- (ww / rskden)[dd == 1]
        rskfactor2 <- (ww / rskden^2)[dd == 1]
        rskdeninv  <- c(0, cumsum(rskfactor))  # pad with 0 so that we can index
        rskdeninv2 <- c(0, cumsum(rskfactor2))

        # this says how many of the risk sets each observation is in; 0 is none
        # (however, if start time is not zero, then we could be including an
        # observation in too many risk sets: we will remove that later.)
        rskcount <- cumsum(dd)
        rskterm  <- rskdeninv[rskcount+1]
        rskterm2 <- rskdeninv2[rskcount+1]
        current_sum <- 0; current_sum2 <- 0
        death_idx <- 1; start_idx <- 1
        while (death_idx <= ndeaths && start_idx <= nobs) {
            if (start_time_start[start_idx] < death_time[death_idx]) {
                # current observation belongs in risk set ending in death time,
                # so we should remove the current cumulative sum and consider
                # the next observation
                stop_idx <- ss_match[start_idx]  # match(start_o[start_idx], stop_o)
                rskterm[stop_idx]  <- rskterm[stop_idx]  - current_sum
                rskterm2[stop_idx] <- rskterm2[stop_idx] - current_sum2
                start_idx <- start_idx + 1
            } else {
                # current observation doesn't belong in risk set ending in death
                # time, so we should add the rskfactor associated with this
                # death time to current_sum and check if the term assoc. with
                # the death time after it should also be added
                current_sum  <- current_sum  + rskfactor[death_idx]
                current_sum2 <- current_sum2 + rskfactor2[death_idx]
                death_idx <- death_idx + 1
            }
        }
        grad <- w * (d - exp_eta * rskterm)
        grad[stop_o] <- grad

        # if diag.hessian = TRUE, return the diagonal of the hessian too
        if (diag.hessian) {
            w_exp_eta <- w * exp_eta
            diag_hessian <- w_exp_eta^2 * rskterm2 - w_exp_eta * rskterm
            diag_hessian[stop_o] <- diag_hessian
            attr(grad, "diag_hessian") <- diag_hessian
        }
        return(grad)
    } else {
        # more than one strata provided: compute strata-level values and
        # concatenate
        overall_grad <- rep(NA, nobs)
        if (diag.hessian) overall_diag_hessian <- rep(NA, nobs)
        for (i in unique(strata)) {
            ii <- which(strata == i)
            strata_res <- coxgrad3(eta[ii], y[ii, , drop = FALSE], w[ii],
                                   std.weights = FALSE, diag.hessian = diag.hessian)
            overall_grad[ii] <- strata_res
            if (diag.hessian) {
                overall_diag_hessian[ii] <- attr(strata_res, "diag_hessian")
            }
        }
        if (diag.hessian) {
            attr(overall_grad, "diag_hessian") <- overall_diag_hessian
        }
        return(overall_grad)
    }
}

#' Helper function for Cox deviance and gradient
#'
#' Helps to find ties in death times of data.
#'
#' @param x Sorted vector of death times.
#' @param index Vector of indices for the death times.
#'
#' @return A list with two arguments.
#' \item{index_first}{A vector of indices for the first observation at each
#' death time as they appear in the sorted list.}
#' \item{index_ties}{If there are no ties at all, this is NULL. If not, this is
#' a list with length equal to the number of unique times with ties. For each
#' time with ties, index_ties gives the indices of the observations with a
#' death at that time.}
#'
#' @examples
#' # Example with no ties
#' glmnet:::fid(c(1, 4, 5, 6), 1:5)
#'
#' # Example with ties
#' glmnet:::fid(c(1, 1, 1, 2, 3, 3, 4, 4, 4), 1:9)
fid <- function(x,index) {
    idup=duplicated(x)
    if(!any(idup)) list(index_first=index,index_ties=NULL)
    else {
        ndup=!idup
        xu=x[ndup]# first death times
        index_first=index[ndup]
        ities=match(x,xu)
        index_ties=split(index,ities)
        nties=sapply(index_ties,length)
        list(index_first=index_first,index_ties=index_ties[nties>1])
    }
}

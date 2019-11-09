#' compute deviance for cox model output
#'
#' Given a fit or coefficients, compute the deciance (-2 log partial likelihood) for
#' right-censored survival data
#'
#' \code{coxnet.deviance} computes the deviance for a single prediction, or a matrix of predictions
#'
#' @aliases coxnet.deviance
#' @param pred matrix of predictions
#' @param y a survival response matrix, as produced by \code{Surv}
#' @param x optional \code{x} matrix, if \code{pred} is \code{NULL}
#' @param offset optional offset
#' @param weights optional observation weights
#' @param beta optional coefficient vector/matrix, supplied if \code{pred=NULL}
#' @return a single or vector of deviances
#' @author Trevor Hastie\cr Maintainer: Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{coxgrad}
#' @keywords Cox model
#'
#' @export coxnet.deviance
coxnet.deviance <-
    function (pred = NULL, y, x = 0, offset = NULL, weights = NULL,
              beta = NULL)
{
    storage.mode(x) = "double"
    y = response.coxnet(y)
    ty = y$time
    tevent = y$event
    ty = ty + (1 - tevent) * 100 * .Machine$double.eps
    nobs = as.integer(length(ty))
    nvars = as.integer(ncol(x))
    nvec=1
    if (is.null(weights))
        weights = rep(1, nobs)
    else {
        weights=nobs*weights/sum(weights)
        weights = as.double(weights)
    }
### Compute saturated loglikelihood
    wd=weights[tevent==1]
    tyd=ty[tevent==1]
    if(any(duplicated(tyd))){
        wd=tapply(wd,tyd,sum)
    }
    wd=wd[wd>0]
    lsat=-sum(wd*log(wd))
####
    if (is.null(offset))
        offset = rep(0, nobs)
    else offset=as.double(offset)
    if (is.null(beta)) {
        beta = double(0)
        nvars = as.integer(0)
    }
    else {
        beta = as.matrix(beta)
        nvec = ncol(beta)
    }
    if(!is.null(pred)){
        # trick to get a set of deviances based on predictions"
        x=as.matrix(pred)
        nvec=ncol(x)
        storage.mode(x)="double"
        beta=diag(nvec)
        nvars=as.integer(nvec)
        storage.mode(beta)="double"
    }
    nvec=as.integer(nvec)

    fit = .Fortran("loglike", nobs, nvars, x, ty, tevent, offset,
                   weights, nvec, beta, flog = double(nvec), jerr = integer(1),
                   PACKAGE = "glmnet")
    if (fit$jerr != 0) {
        errmsg = jerr(fit$jerr, maxit = 0, pmax = 0, family = "cox")
        if (errmsg$fatal)
            stop(errmsg$msg, call. = FALSE)
        else warning(errmsg$msg, call. = FALSE)
    }
    2 *(lsat-fit$flog)
}

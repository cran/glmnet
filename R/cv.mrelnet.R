cv.mrelnet <-
    function (outlist, lambda, x, y, weights, offset, foldid, type.measure,
              grouped, keep = FALSE, alignment="lambda")
{
    ndim = dim(y)
    nc = ndim[2]
    nobs = ndim[1]
    if (!is.null(offset))
        y = y - drop(offset)
    predmat = array(NA, c(nobs, nc, length(lambda)))
    nfolds = max(foldid)
    nlams = double(nfolds)
    nlambda=length(lambda)
    for (i in seq(nfolds)) {
        which = foldid == i
        fitobj = outlist[[i]]
        fitobj$offset = FALSE
        preds = switch(alignment,
                       fraction = predict(fitobj, x[which, , drop = FALSE]),
                                          lambda = predict(fitobj, x[which, , drop = FALSE], s = lambda)
                                          )
        nlami = min(dim(preds)[3],nlambda)
        predmat[which, , seq(nlami)] = preds[,,seq(nlami)]
        if(nlami<nlambda)predmat[which,,seq(from=nlami,to=nlambda)]=preds[,,nlami]
        nlams[i] = nlambda
    }
    N = nobs - apply(is.na(predmat[, 1, ]), 2, sum)
    bigY = array(y, dim(predmat))
    cvraw = switch(type.measure, mse = apply((bigY - predmat)^2,
                                             c(1, 3), sum), deviance = apply((bigY - predmat)^2, c(1,
                                                                                                   3), sum), mae = apply(abs(bigY - predmat), c(1, 3), sum))
    if ((nobs/nfolds < 3) && grouped) {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",
                call. = FALSE)
        grouped = FALSE
    }
    if (grouped) {
        cvob = cvcompute(cvraw, weights, foldid, nlams)
        cvraw = cvob$cvraw
        weights = cvob$weights
        N = cvob$N
    }
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                      w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, type.measure = type.measure)
    if (keep)
        out$fit.preval = predmat
    out
}

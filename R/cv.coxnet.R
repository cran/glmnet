cv.coxnet <- function(predmat,y,type.measure,weights,foldid,grouped){
    ## Note, all the work got done in buildPredmat.coxnetlist for deviance; a special case
devtrue=type.measure=="deviance"
    nfolds = max(foldid)

    if ((length(weights)/nfolds < 10) && !grouped) {
        warning("Option grouped=TRUE enforced for cv.coxnet, since < 10 observations per fold",
                call. = FALSE)
    }

    if(devtrue){
        cvraw=attr(predmat,"cvraw")
        status = y[, "status"]
        N = nfolds - apply(is.na(cvraw), 2, sum)
        weights = as.vector(tapply(weights * status, foldid, sum))
        cvraw = cvraw/weights
    }
    else
    {
        nlambda=ncol(predmat)
        nlams=rep(nlambda,nfolds)
        cvraw = matrix(NA, nfolds, nlambda)
        good = matrix(0, nfolds, nlambda)
        for (i in seq(nfolds)) {
            good[i, seq(nlams[i])] = 1
            which = foldid == i
            for (j in seq(nlams[i])) {
                cvraw[i, j] = Cindex(predmat[which,j],y[which, ], weights[which])
            }
        }
        N = apply(good, 2, sum)
        weights = tapply(weights, foldid, sum)
        }
    list(cvraw=cvraw,weights=weights,N=N,type.measure=type.measure,grouped=FALSE)

}

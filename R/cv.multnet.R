cv.multnet <-function(predmat,y,type.measure,weights,foldid,grouped){
    prob_min = 1e-05
    prob_max = 1 - prob_min
    nc = dim(y)
    if (is.null(nc)) {
        y = as.factor(y)
        ntab = table(y)
        nc = as.integer(length(ntab))
        y = diag(nc)[as.numeric(y), ,drop=FALSE]
    }
    else nc = nc[2]

    ywt = apply(y, 1, sum)
    y = y/ywt
    weights = weights * ywt
    N = nrow(y) - apply(is.na(predmat[, 1, ,drop=FALSE]), 2, sum)## dimensions could be lost if third dim=1
    bigY = array(y, dim(predmat))
    predmat=exp(predmat)
    predtot=apply(predmat,c(1,3),sum)
    for(i in 1:nc)predmat[,i,]=predmat[,i,]/predtot
    cvraw = switch(type.measure,
                   mse = apply((bigY - predmat)^2,c(1, 3), sum),
                   mae = apply(abs(bigY - predmat), c(1,3), sum),
                   deviance = {                                                                                                              predmat = pmin(pmax(predmat, prob_min), prob_max)
                         lp = bigY * log(predmat)
                         ly = bigY * log(bigY)
                         ly[bigY == 0] = 0
                         apply(2 * (ly - lp), c(1, 3), sum)
                   },
                   class = {
                       classid = as.numeric(apply(predmat, 3, glmnet_softmax,ignore_labels=TRUE))
                       yperm = matrix(aperm(bigY, c(1, 3, 2)), ncol = nc)
                       matrix(1 - yperm[cbind(seq(classid), classid)], ncol = ncol(predtot))
                   }
                   )
    list(cvraw=cvraw,weights=weights,N=N,type.measure=type.measure,grouped=grouped)
}



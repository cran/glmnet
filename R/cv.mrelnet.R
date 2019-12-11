cv.mrelnet <- function(predmat,y,type.measure,weights,foldid,grouped){
    ndim = dim(y)
    nc = ndim[2]
    nobs = ndim[1]
    N = nobs - apply(is.na(predmat[, 1,,drop=FALSE ]), 2, sum)# dimensions could be lost if third dim=1
    bigY = array(y, dim(predmat))
    cvraw = switch(type.measure,
                   mse = apply((bigY - predmat)^2,c(1, 3), sum),
                   deviance = apply((bigY - predmat)^2, c(1,3), sum),
                   mae = apply(abs(bigY - predmat), c(1, 3), sum)
                   )
    list(cvraw=cvraw,weights=weights,N=N,type.measure=type.measure,grouped=grouped)
}

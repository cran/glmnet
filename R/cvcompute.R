cvcompute <-
    function (cvstuff, foldid, nlams)
{
    weights=cvstuff$weights
    mat=cvstuff$cvraw
    wisum = tapply(weights, foldid, sum)
    nfolds = max(foldid)
    outmat = matrix(NA, nfolds, ncol(mat))
    good = matrix(0, nfolds, ncol(mat))
    mat[is.infinite(mat)] = NA
    for (i in seq(nfolds)) {
        mati = mat[foldid == i, , drop = FALSE]
        wi = weights[foldid == i]
        outmat[i, ] = apply(mati, 2, weighted.mean, w = wi, na.rm = TRUE)
        good[i, seq(nlams[i])] = 1
    }
    N = apply(good, 2, sum)
    list(cvraw = outmat, weights = wisum, N = N, type.measure=cvstuff$type.measure)
}

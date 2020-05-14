#' @method buildPredmat array
#' @export
buildPredmat.array=function(outlist, lambda, x, offset, foldid, alignment,...){
    if (!is.null(offset)) is.offset = TRUE
    else is.offset = FALSE
    nc=dim(outlist[[1]]$a0)[1]
    nobs=nrow(x)
    predmat = array(NA, c(nobs, nc, length(lambda)))
    nfolds = max(foldid)
    nlams = double(nfolds)
    nlambda = length(lambda)
    for (i in seq(nfolds)) {
        which = foldid == i
        fitobj = outlist[[i]]
        if (is.offset)
            off_sub = offset[which, , drop = FALSE]
        preds = switch(alignment, fraction = predict(fitobj,
                       x[which, , drop = FALSE], newoffset = off_sub,...),
                       lambda = predict(fitobj, x[which, , drop = FALSE],
                       s = lambda, newoffset = off_sub,...))
        nlami = min(dim(preds)[3], nlambda)
        predmat[which, , seq(nlami)] = preds[, , seq(nlami)]
        if (nlami < nlambda)
            predmat[which, , seq(from = nlami, to = nlambda)] = preds[,, nlami]
    }
### fix up dimnames
    rn=rownames(x)
    sn=paste("s",seq(0,length=nlambda),sep="")
    cn=dimnames(preds)[[2]]
    dimnames(predmat)=list(rn,cn,sn)
    predmat
}

#' @method buildPredmat multnetlist
#' @export
buildPredmat.multnetlist=function(outlist, lambda, x, offset, foldid, alignment,...){
    mat=buildPredmat.array(outlist, lambda, x, offset, foldid, alignment,...)
    attr(mat,"classnames")=outlist[[1]]$classnames
    mat
}

#' @method buildPredmat mrelnetlist
#' @export
buildPredmat.mrelnetlist=function(outlist, lambda, x, offset, foldid, alignment,...){
    buildPredmat.array(outlist, lambda, x, offset, foldid, alignment,...)
}

#' @method buildPredmat lognetlist
#' @export
buildPredmat.lognetlist=function(outlist, lambda, x, offset, foldid, alignment,...){
    mat=buildPredmat.default(outlist, lambda, x, offset, foldid, alignment,...)
    attr(mat,"classnames")=outlist[[1]]$classnames
    mat
}


#' @method buildPredmat glmnetfitlist
#' @export
buildPredmat.glmnetfitlist=function(outlist, lambda, x, offset, foldid, alignment,family,...){
    mat=buildPredmat.default(outlist, lambda, x, offset, foldid, alignment,...)
    attr(mat,"family")=family
    mat
}

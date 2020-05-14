#' @export
buildPredmat.default=function(outlist, lambda, x, offset, foldid, alignment,...){
    if (!is.null(offset)) {
        is.offset = TRUE
        offset = drop(offset)
    }
    else is.offset = FALSE
    predmat = matrix(NA, nrow(x), length(lambda))
    nfolds = max(foldid)
    nlams = double(nfolds)
    nlambda=length(lambda)
    for (i in seq(nfolds)) {
        which = foldid == i
        fitobj = outlist[[i]]
        if (is.offset)
            off_sub = offset[which]
        preds = switch(alignment,
                   fraction=predict(fitobj, x[which, , drop = FALSE], newoffset = off_sub,...),
                   lambda=predict(fitobj, x[which, , drop = FALSE], s=lambda, newoffset = off_sub,...)
                   )
    nlami = min(ncol(preds),nlambda)
    predmat[which, seq(nlami)] = preds[,seq(nlami)]
    if(nlami<nlambda)predmat[which,seq(from=nlami,to=nlambda)]=preds[,nlami]
    }
    rn=rownames(x)
    sn=paste("s",seq(0,length=nlambda),sep="")
    dimnames(predmat)=list(rn,sn)

    predmat
    }

#' @export
buildPredmat=function(outlist, lambda, x, offset, foldid, alignment,...)
    UseMethod("buildPredmat")


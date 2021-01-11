#' @method buildPredmat coxnetlist
#' @export
buildPredmat.coxnetlist <-
    function(outlist, lambda, x, offset, foldid, alignment,y,weights,grouped,type.measure="deviance",...){
        nfolds = max(foldid)
        if ((length(weights)/nfolds < 10) && !grouped) grouped = TRUE
        devtrue=type.measure=="deviance"
        cvraw = if(devtrue) matrix(NA, nfolds, length(lambda)) else NULL
        nlambda = length(lambda)
        predmat=matrix(NA, nrow(x), length(lambda))
        rn=rownames(x)
        sn=paste("s",seq(0,length=nlambda),sep="")
        dimnames(predmat)=list(rn,sn)
    predmat

    for (i in seq(nfolds)) {
        which = foldid == i
        fitobj = outlist[[i]]
        coefmat = switch(alignment,
                         fraction = predict(fitobj,type = "coefficients",...),
                         lambda = predict(fitobj, type = "coefficients",s = lambda,...)
                         )
        nlami = min(ncol(coefmat), nlambda)
        if(devtrue){
            if (grouped) {
                plfull = coxnet.deviance(x = x, y = y, offset = offset,
                                         weights = weights, beta = coefmat)
                plminusk = coxnet.deviance(x = x[!which, ], y = y[!which, ], 
                                           offset = offset[!which],
                                           weights = weights[!which], 
                                           beta = coefmat)
                cvraw[i, seq(nlami)] = (plfull - plminusk)[seq(nlami)]
            }
            else {
                plk = coxnet.deviance(x = x[which, ], y = y[which, ], 
                                      offset = offset[which], 
                                      weights = weights[which], beta = coefmat)
                cvraw[i, seq(nlami)] = plk[seq(nlami)]
            }
        }
        predmat[which, seq(nlami)] = as.matrix(x[which,]%*% coefmat)
        if (nlami < nlambda){
            if(devtrue)cvraw[i, seq(from = nlami, to = nlambda)] = cvraw[i, nlami]
            predmat[which,seq(from=nlami,to=nlambda)]=predmat[which,nlami]
        }
    }
        if(devtrue) attr(predmat,"cvraw")=cvraw
        predmat
    }

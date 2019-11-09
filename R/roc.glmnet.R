#' Compute the roc curve(s) for a binimial glmnet model.
#'
#' @rdname assess.glmnet
#' @export roc.glmnet
roc.glmnet=function(object,newx=NULL,newy,...){
    ### object must be either a glmnet or cv.glmnet object, or else a matrix or a vector of predictions of a glmnet model
    oclass=grep("glmnet",class(object))
    if(length(oclass)){# either a glmnet or cv.glmnet
        if(family(object)!="binomial")stop("roc available only for binomial family")
        predmat=predict(object,newx=newx,...)
    }
    else{
        predmat=object
        if(is.vector(predmat))predmat=as.matrix(predmat)
        }

roc.raw=function(prediction,y){
    op=order(prediction,decreasing=TRUE)
    preds=prediction[op]
    y=y[op]
    noty=1-y
    if(any(duplicated(preds))){
        y=rev(tapply(y,preds,sum))
        noty=rev(tapply(noty,preds,sum))
    }
   data.frame( FPR=cumsum(noty)/sum(noty),TPR=cumsum(y)/sum(y))
}

### newy should be binary, but could be a two-class matrix!
      nc = dim(newy)
    if (is.null(nc)) {
        newy = as.factor(newy)
        ntab = table(newy)
        nc = length(ntab)
        newy = diag(nc)[as.numeric(newy), ]
    }
    else {
         nc = nc[2]
    }
    if(nc!=2)stop("newy should be binary, or a factor with two levels")
    newy=newy[,2]
    n=nrow(predmat)
    if(n!=length(newy))stop("mismatch in size of two arguments")
    out=apply(predmat,2,roc.raw,y=newy)
    if(length(out)==1)out=out[[1]]
    out
}




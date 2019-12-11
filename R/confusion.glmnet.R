#' Compute a confusion matrix (or many) for a binomial or multinomial glmnet model
#'
#' @rdname assess.glmnet
#' @export confusion.glmnet
confusion.glmnet=function(object,newx=NULL,newy,family=c("binomial","multinomial"),...){
### object must be either a glmnet or cv.glmnet object fit with family binomial or multinomial
### or else a matrix/array of predictions of a glmnet model of these families
###    (the last dimension can be 1)
###    (not on the mean scale, but natural parameter scale)
    oclass=grep("glmnet",class(object))
    if(length(oclass)){# either a glmnet or cv.glmnet
        if(!match(family(object),c("binomial","multinomial"),FALSE))
            stop("confusion available only for binomial or multinomial family")
        classmat=predict(object,newx=newx,type="class",...)
    }
    else{
        fam=match.arg(family)
        predmat=object
        if(fam=="binomial"){
            classnames=attr(predmat,"classnames")#if it was created by cv with keep
            if(is.null(classnames))classnames=c("1","2")
            if(is.vector(predmat))predmat=as.matrix(predmat)
            cnum=ifelse(predmat>0,2,1)
            classmat=classnames[cnum]
            if(is.matrix(predmat))classmat=array(classmat,dim(predmat),dimnames(predmat))
        }
        else{
            classmat=apply(predmat, 3, glmnet_softmax)
        }
    }

### newy can be a matrix or a factor
      nc = dim(newy)
    if (is.null(nc)) {
        newy = as.factor(newy)
    }
    else {
        nc = nc[2]
        yc=newy%*%(1:nc)
        cn=colnames(newy)
        if(!is.null(cn))
            newy=factor(yc,labels=cn)
        else newy=factor(yc)
    }
    ctable=function(...){
        tab=table(...)
        class(tab)=c("confusion.table",class(tab))
        tab
        }
    if(ncol(classmat)>1)
        ## convert to a dataframe, to prevent apply simplifying
        lapply(data.frame(classmat,stringsAsFactors=FALSE,check.names=FALSE),
               function(x,y)ctable(Predicted=x,True=y),y=newy)
    else(ctable(Predicted=classmat,True=newy))
}




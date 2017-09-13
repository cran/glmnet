check_dots<-
    function(object,...){
        thiscall=object$call
        ncall=names(thiscall)[-1]
        needarg=c("x","y","weights","offset","penalty.factor","lower.limits","upper.limits")
        w=match(ncall,needarg,0)
        needarg=needarg[w]
        nargs=names(list(...))
        w=match(needarg,nargs,0)>0
        if(!all(w)){
            margs=needarg[!w]
                stop(paste("used coef.glmnet() or predict.glmnet() with `exact=TRUE` so must in addition supply original argument(s) ",paste(margs,collapse=" and "), " in order to safely rerun glmnet"),call.=FALSE)
        }
        invisible()
    }



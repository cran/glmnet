check_dots<-
    function(object,...,need=c("x","y","weights","offset","penalty.factor","lower.limits","upper.limits"),error="used coef.glmnet() or predict.glmnet() with `exact=TRUE`"){
        if(is.null(need))return(invisible())
        thiscall=object$call
        ncall=names(thiscall)[-1]
        w=match(ncall,need,0)
        need=need[w]
        nargs=names(list(...))
        w=match(need,nargs,0)>0
        if(!all(w)){
            margs=need[!w]
                stop(paste(error,"so must in addition supply original argument(s) ",paste(margs,collapse=" and "), " in order to safely rerun glmnet"),call.=FALSE)
        }
        invisible()
    }



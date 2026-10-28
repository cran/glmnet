#' @method predict relaxed
#' @param gamma Value(s) of \code{gamma} at which predictions are required,
#' for "relaxed" objects. If more than one value is given, a list of
#' predictions is returned, one per \code{gamma}, named by the rounded
#' \code{gamma} values.
#' @rdname predict.glmnet
#' @export
#' @export predict.relaxed
predict.relaxed=function (object, newx, s = NULL, gamma=1,type = c("link", "response",
    "coefficients", "nonzero", "class"), exact = FALSE, newoffset,
    ...) {
    exact=FALSE # we cannot accommodate exact here
    type=match.arg(type)
    gamma=checkgamma.relax(gamma)
    if(length(gamma)>1){
        ng=length(gamma)
        outlist=vector("list", ng)
        for( i in 1:ng)outlist[[i]]=predict(object, newx, s, gamma=gamma[i], type=type,
                                            exact=exact, newoffset=newoffset, ...)
        names(outlist)=format(round(gamma,2))
        return(outlist)
    }
    if(gamma==1)return(NextMethod("predict"))
    predict(blend.relaxed(object,gamma),newx, s,type, exact, newoffset,  ...)
}

checkgamma.relax=function(gamma){
    if(any(wh<-gamma<0)){
        warning("negative gamma values ignored")
        gamma=gamma[!wh]
    }
    if(any(wh<-gamma>1)){
        warning("gamma values larger than 1 ignored")
        gamma=gamma[!wh]
    }
    if(!length(gamma))stop("no valid values of gamma")
    gamma
}


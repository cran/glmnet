check_dots <- function(object, ..., 
                       need = c("x", "y", "weights", "offset", "penalty.factor",
                                "lower.limits", "upper.limits"), 
                       error_start = "used coef.glmnet() or predict.glmnet() with `exact=TRUE`",
                       error_end = " in order to safely rerun glmnet",
                       prefix = NULL) {
    if (is.null(need)) return(invisible())
    
    # extract the function options we need from the object's call
    thiscall = object$call
    ncall = names(thiscall)[-1]
    w = match(ncall, need, 0)
    need = need[w]
    if (length(need) == 0) return(invisible())
    
    # check that ... indeed has those function options
    if (!is.null(prefix)) need <- paste0(prefix, need)
    nargs = names(list(...))
    w = match(need, nargs, 0) > 0
    if(!all(w)) {
        margs = need[!w]
        stop(paste(error_start, 
                   "so must in addition supply original argument(s) ",
                   paste(margs,collapse=" and "), 
                   error_end), call.=FALSE)
    }
    invisible()
}


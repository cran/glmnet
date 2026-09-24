check.penalty.factor <- function(penalty.factor,nvars){
    if(length(penalty.factor)!=nvars)
        stop("the length of penalty.factor does not match the number of variables")
    if(any(penalty.factor < 0)){
        warning("values for penalty.factor must be between 0 and Inf; negative values set to zero")
        penalty.factor = pmax(0,penalty.factor)
    }
    penalty.factor
}

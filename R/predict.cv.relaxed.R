#' @method predict cv.relaxed
#' @param gamma Value (single) of 'gamma' at which predictions are to be made
#' @rdname predict.cv.glmnet
#' @export
 predict.cv.relaxed <-
function (object, newx, s = c("lambda.1se", "lambda.min"), gamma=c("gamma.1se","gamma.min"),...)
{
    relaxed=object$relaxed
    if (is.numeric(s))
        lambda = s
    else if (is.character(s)) {
        s = match.arg(s)
        lambda = relaxed[[s]]
        names(lambda)=s
        if(missing(gamma)){
            gamma=switch(s,
                         lambda.1se=relaxed[["gamma.1se"]],
                         lambda.min=relaxed[["gamma.min"]]
                         )
            }
    }
else stop("Invalid form for s")
    if (is.character(gamma)) {
        gamma = match.arg(gamma)
        gamma = relaxed[[gamma]]
    }
    if(!is.numeric(gamma))stop("Invalid form for gamma")

    predict(object$glmnet.fit, newx, s = lambda,gamma=gamma, ...)
}

#' @method coef cv.relaxed
#' @export
coef.cv.relaxed <-
function (object, s = c("lambda.1se", "lambda.min"), gamma=c("gamma.1se","gamma.min"),...)
 {
     relaxed=object$relaxed
    if (is.numeric(s))
        lambda = s
    else if (is.character(s)) {
        s = match.arg(s)
        lambda = relaxed[[s]]
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

    coef(object$glmnet.fit, s = lambda, gamma=gamma, ...)
}

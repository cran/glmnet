#' @method family glmnet
#' @export
family.glmnet=function(object,...){
    families=c(elnet = "gaussian", lognet = "binomial", fishnet = "poisson",
               multnet = "multinomial", coxnet = "cox", mrelnet = "mgaussian")
    cl=class(object)[1]
    families[cl]
}

#' @method family glmnetfit
#' @export
family.glmnetfit=function(object,...) object$family


#' @method family relaxed
#' @export
family.relaxed=function(object,...)family(object$relaxed)

#' @method family cv.glmnet
#' @export
family.cv.glmnet=function(object,...)family(object$glmnet.fit)


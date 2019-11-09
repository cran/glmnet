#' Extract coefficients from a glmnet object
#'
#' @method coef glmnet
#' @rdname predict.glmnet
#' @export
#' @export coef.glmnet
coef.glmnet=function(object,s=NULL,exact=FALSE,...)
  predict(object,s=s,type="coefficients",exact=exact,...)

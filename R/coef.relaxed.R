#' @method coef relaxed
#' @export
#' @export coef.relaxed
coef.relaxed=function(object,s=NULL, gamma=1,...)
    predict(object, s=s, gamma= gamma, type="coefficients",...)

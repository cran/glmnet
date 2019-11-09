#' @method predict elnet
#' @export
predict.elnet=function(object,newx,s=NULL,type=c("link","response","coefficients","nonzero"),exact=FALSE,newoffset,...){
NextMethod("predict")
 }

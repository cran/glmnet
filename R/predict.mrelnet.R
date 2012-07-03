predict.mrelnet=function(object, newx, s = NULL, type = c("link", "response", 
                                      "coefficients", "nonzero"), exact = FALSE, offset, ...)
  {
  type=match.arg(type)
  if(type=="response")type="link"
  object$grouped=TRUE
   predict.multnet(object,newx,s,type,exact,offset,...)
}

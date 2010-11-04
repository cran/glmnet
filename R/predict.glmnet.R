predict.glmnet=function(object,newx,s=NULL,type=c("link","response","coefficients","class","nonzero"),exact=FALSE,offset,...){
  type=match.arg(type)
  if(missing(newx)){
    if(!match(type,c("coefficients","nonzero"),FALSE))stop("You need to supply a value for 'newx'")
     }
NextMethod("predict")
  }

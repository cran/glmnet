predict.cv.glmnet=function(object,newx,s=c("lambda.1se","lambda.min"),...){
  if(is.numeric(s))lambda=s
  else 
    if(is.character(s)){
      s=match.arg(s)
      lambda=object[[s]]
    }
    else stop("Invalid form for s")
  predict(object$glmnet.fit,newx,s=lambda,...)
}

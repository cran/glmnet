plot.glmnet=function(x, xvar=c("norm","lambda","dev"),label=FALSE,...){
  xvar=match.arg(xvar)
 plotCoef(x$beta,lambda=x$lambda,df=x$df,dev=x$dev.ratio,label=label,xvar=xvar,...)
}

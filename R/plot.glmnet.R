plot.glmnet=function(x, xvar=c("norm","lambda","dev"),label=FALSE,...){
  ocl=class(x)[[2]]
  xvar=match.arg(xvar)
  if(ocl=="multnet"){
    beta=x$beta
    if(xvar=="norm"){
      cnorm1=function(beta){
        which=nonzeroCoef(beta)
        beta=as.matrix(beta[which,])
        apply(abs(beta),2,sum)
      }
      norm=apply(sapply(x$beta,cnorm1),1,sum)
    } else norm = NULL
    dfmat=x$dfmat
    ncl=nrow(dfmat)
    clnames=rownames(dfmat)
    for( i in seq(ncl)){
      plotCoef(beta[[i]],norm,x$lambda,dfmat[i,],x$dev.ratio,label=label,xvar=xvar,ylab=paste("Coefficients: Class",clnames[i]),...)
    }
  }else plotCoef(x$beta,lambda=x$lambda,df=x$df,dev=x$dev.ratio,label=label,xvar=xvar,...)
}
      
      

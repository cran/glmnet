plotCoef=function(beta,norm,lambda,df,dev,label=FALSE,xvar=c("norm","lambda","dev"),xlab=iname,ylab="Coefficients",...){
  ##beta should be in "dgCMatrix" format
  which=nonzeroCoef(beta)
  beta=as.matrix(beta[which,])
  xvar=match.arg(xvar)
  switch(xvar,
    "norm"={
      index=if(missing(norm))apply(abs(beta),2,sum)else norm
      iname="L1 Norm"
    },
    "lambda"={
      index=log(lambda)
      iname="Log Lambda"
    },
    "dev"= {
      index=dev
      iname="Fraction Deviance Explained"
    }
         )
  dotlist=list(...)
  type=dotlist$type
  if(is.null(type))
    matplot(index,t(beta),lty=1,xlab=xlab,ylab=ylab,type="l",...)
  else matplot(index,t(beta),lty=1,xlab=xlab,ylab=ylab,...)
  atdf=pretty(index)
 prettydf=trunc(approx(x=index,y=df,xout=atdf,rule=2)$y)
 axis(3,at=atdf,label=prettydf,tcl=NA)
 if(label){
   nnz=length(which)
   xpos=max(index)
   pos=4
   if(xvar=="lambda"){
     xpos=min(index)
     pos=2
   }
   xpos=rep(xpos,nnz)
   ypos=beta[,ncol(beta)]
   text(xpos,ypos,paste(which),cex=.5,pos=pos)
 }
         
}
    

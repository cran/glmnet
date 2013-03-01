elnet=function(x,is.sparse,ix,jx,y,weights,offset,type.gaussian=c("covariance","naive"),alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,vnames,maxit){
  maxit=as.integer(maxit)
  weights=as.double(weights)
### compute the null deviance
  ybar=weighted.mean(y,weights)
  nulldev=sum(weights* (y-ybar)^2)
  type.gaussian=match.arg(type.gaussian)

  ka=as.integer(switch(type.gaussian,
    covariance=1,
    naive=2,
    ))


 storage.mode(y)="double"
   if(is.null(offset)){
    offset=y*0
    is.offset=FALSE}
  else{
    storage.mode(offset)="double"
    is.offset=TRUE
  }

fit=if(is.sparse).Fortran("spelnet",
        ka,parm=alpha,nobs,nvars,x,ix,jx,y-offset,weights,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,maxit,
        lmu=integer(1),
        a0=double(nlam),
        ca=double(nx*nlam),
        ia=integer(nx),
        nin=integer(nlam),
        rsq=double(nlam),
        alm=double(nlam),
        nlp=integer(1),
        jerr=integer(1),PACKAGE="glmnet"
        )
else .Fortran("elnet",
          ka,parm=alpha,nobs,nvars,as.double(x),y-offset,weights,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,maxit,
          lmu=integer(1),
          a0=double(nlam),
          ca=double(nx*nlam),
          ia=integer(nx),
          nin=integer(nlam),
          rsq=double(nlam),
          alm=double(nlam),
          nlp=integer(1),
          jerr=integer(1),PACKAGE="glmnet"
          )
if(fit$jerr!=0){
  errmsg=jerr(fit$jerr,maxit,pmax=nx,family="gaussian")
  if(errmsg$fatal)stop(errmsg$msg,call.=FALSE)
  else warning(errmsg$msg,call.=FALSE)
}
  outlist=getcoef(fit,nvars,nx,vnames)
  dev=fit$rsq[seq(fit$lmu)]
  outlist=c(outlist,list(dev.ratio=dev,nulldev=nulldev,npasses=fit$nlp,jerr=fit$jerr,offset=is.offset))
  class(outlist)="elnet"
  outlist
}

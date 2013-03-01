fishnet=function(x,is.sparse,ix,jx,y,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,vnames,maxit){
  if(any(y<0))stop("negative responses encountered;  not permitted for Poisson family")
  maxit=as.integer(maxit)
  weights=as.double(weights)
  storage.mode(y)="double"
   if(is.null(offset)){
    offset=y*0 #keeps the shape of y
    is.offset=FALSE}
  else{
    storage.mode(offset)="double"
    
    is.offset=TRUE
  }
fit=if(is.sparse).Fortran("spfishnet",
  parm=alpha,nobs,nvars,x,ix,jx,y,offset,weights,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,maxit,
  lmu=integer(1),
  a0=double(nlam),
  ca=double(nx*nlam),
  ia=integer(nx),
  nin=integer(nlam),
  nulldev=double(1),
  dev=double(nlam),
  alm=double(nlam),
  nlp=integer(1),
  jerr=integer(1),PACKAGE="glmnet"
  )
else .Fortran("fishnet",
              parm=alpha,nobs,nvars,as.double(x),y,offset,weights,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,maxit,
              lmu=integer(1),
              a0=double(nlam),
              ca=double(nx*nlam),
              ia=integer(nx),
              nin=integer(nlam),
              nulldev=double(1),
              dev=double(nlam),
              alm=double(nlam),
              nlp=integer(1),
              jerr=integer(1),PACKAGE="glmnet"
              )
if(fit$jerr!=0){
  errmsg=jerr(fit$jerr,maxit,pmax=nx,family="poisson")
  if(errmsg$fatal)stop(errmsg$msg,call.=FALSE)
  else warning(errmsg$msg,call.=FALSE)
}
  outlist=getcoef(fit,nvars,nx,vnames)
  dev=fit$dev[seq(fit$lmu)]
outlist=c(outlist,list(dev.ratio=dev,nulldev=fit$nulldev,npasses=fit$nlp,jerr=fit$jerr,offset=is.offset))
class(outlist)="fishnet"
outlist
}

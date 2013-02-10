coxnet=function(x,is.sparse,ix,jx,y,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,vnames,maxit){
  if(!is.matrix(y)||!all(match(c("time","status"),dimnames(y)[[2]],0)))stop("Cox model requires a matrix with columns 'time' (>0) and 'status'  (binary) as a response; a 'Surv' object suffices",call.=FALSE)
  ty=as.double(y[,"time"])
  tevent=as.double(y[,"status"])
  if(any(ty<=0))stop("negative event times encountered;  not permitted for Cox family")
  maxit=as.integer(maxit)
  weights=as.double(weights)
     if(is.null(offset)){
    offset=ty*0 
    is.offset=FALSE}
  else{
    storage.mode(offset)="double"
    is.offset=TRUE
  }
  fit=if(is.sparse)stop("Cox model mot implemented for sparse x in glmnet",call.=FALSE)
  else .Fortran("coxnet",
              parm=alpha,nobs,nvars,as.double(x),ty,tevent,offset,weights,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,
              maxit,isd,# need to get JHF to reverse these
              lmu=integer(1),
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
  errmsg=jerr(fit$jerr,maxit,pmax=nx,family="cox")
  if(errmsg$fatal)stop(errmsg$msg,call.=FALSE)
  else warning(errmsg$msg,call.=FALSE)
}
  outlist=getcoef(fit,nvars,nx,vnames)
  dev=fit$dev[seq(fit$lmu)]
outlist=c(outlist,list(dev.ratio=dev,nulldev=fit$nulldev,npasses=fit$nlp,jerr=fit$jerr,offset=is.offset))
class(outlist)="coxnet"
outlist
}

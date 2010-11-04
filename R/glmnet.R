glmnet=function(x,y,family=c("gaussian","binomial","poisson","multinomial","cox"),weights,offset=NULL,alpha=1.0,nlambda=100,lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4),lambda=NULL,standardize=TRUE,thresh=1e-4,dfmax=nvars+1,pmax=min(dfmax*1.2,nvars),exclude,penalty.factor=rep(1,nvars),maxit=100,HessianExact=FALSE,type=c("covariance","naive")){

### Prepare all the generic arguments, then hand off to family functions
  family=match.arg(family)
  this.call=match.call()
  nlam=as.integer(nlambda)
  y=drop(y) # we dont like matrix responses unless we need them
  np=dim(x)
  nobs=as.integer(np[1])
  if(missing(weights))weights=rep(1,nobs)
  nvars=as.integer(np[2])
  vnames=colnames(x)
  if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")
  ne=as.integer(dfmax)
  nx=as.integer(pmax)
  if(!missing(exclude)){
    jd=match(exclude,seq(nvars),0)
    if(!all(jd>0))stop("Some excluded variables out of range")
    jd=as.integer(c(length(jd),jd))
  }else jd=as.integer(0)
  vp=as.double(penalty.factor)
  isd=as.integer(standardize)
  thresh=as.double(thresh)
  if(is.null(lambda)){
     if(lambda.min.ratio>=1)stop("lambda.min.ratio should be less than 1")
    flmin=as.double(lambda.min.ratio)
    ulam=double(1)
  }
  else{
     flmin=as.double(1)    
    if(any(lambda<0))stop("lambdas should be non-negative")
    ulam=as.double(rev(sort(lambda)))
    nlam=as.integer(length(lambda))
  }
  is.sparse=FALSE
  ix=jx=NULL
  if(class(x)=="dgCMatrix"){##Sparse case
    is.sparse=TRUE
    ix=as.integer(x@p+1)
    jx=as.integer(x@i+1)
    x=as.double(x@x)
  }

  fit=switch(family,
    "gaussian"=elnet(x,is.sparse,ix,jx,y,weights,offset,type,alpha,nobs,nvars,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,vnames),
    "poisson"=fishnet(x,is.sparse,ix,jx,y,weights,offset,type,alpha,nobs,nvars,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,vnames,maxit),
    "binomial"=lognet(x,is.sparse,ix,jx,y,weights,offset,type,alpha,nobs,nvars,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,vnames,maxit,!HessianExact,family),
    "multinomial"=lognet(x,is.sparse,ix,jx,y,weights,offset,type,alpha,nobs,nvars,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,vnames,maxit,!HessianExact,family),
    "cox"=coxnet(x,is.sparse,ix,jx,y,weights,offset,type,alpha,nobs,nvars,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,vnames,maxit)
    )
    
  if(is.null(lambda))fit$lambda=fix.lam(fit$lambda)##first lambda is infinity; changed to entry point
fit$call=this.call
  class(fit)=c("glmnet",class(fit))
  fit
}

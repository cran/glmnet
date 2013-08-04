glmnet=function(x,y,family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),weights,offset=NULL,alpha=1.0,nlambda=100,lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4),lambda=NULL,standardize=TRUE,intercept=TRUE,thresh=1e-7,dfmax=nvars+1,pmax=min(dfmax*2+20,nvars),exclude,penalty.factor=rep(1,nvars),lower.limits=-Inf,upper.limits=Inf,maxit=100000,type.gaussian=ifelse(nvars<500,"covariance","naive"),type.logistic=c("Newton","modified.Newton"),standardize.response=FALSE,type.multinomial=c("ungrouped","grouped")){

### Prepare all the generic arguments, then hand off to family functions
  family=match.arg(family)
  if(alpha>1){
    warning("alpha >1; set to 1")
    alpha=1
  }
  if(alpha<0){
    warning("alpha<0; set to 0")
    alpha=0
  }
  alpha=as.double(alpha)

  this.call=match.call()
  nlam=as.integer(nlambda)
  y=drop(y) # we dont like matrix responses unless we need them
  np=dim(x)
   ###check dims
  if(is.null(np)|(np[2]<=1))stop("x should be a matrix with 2 or more columns")
 nobs=as.integer(np[1])
  if(missing(weights))weights=rep(1,nobs)
  else if(length(weights)!=nobs)stop(paste("number of elements in weights (",length(weights),") not equal to the number of rows of x (",nobs,")",sep=""))
  nvars=as.integer(np[2])
  dimy=dim(y)
  nrowy=ifelse(is.null(dimy),length(y),dimy[1])
    if(nrowy!=nobs)stop(paste("number of observations in y (",nrowy,") not equal to the number of rows of x (",nobs,")",sep=""))
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
 ###check on limits
  internal.parms=glmnet.control()
  if(any(lower.limits>0)){stop("Lower limits should be non-positive")}
  if(any(upper.limits<0)){stop("Upper limits should be non-negative")}
  lower.limits[lower.limits==-Inf]=-internal.parms$big
  upper.limits[upper.limits==Inf]=internal.parms$big
  if(length(lower.limits)<nvars){
    if(length(lower.limits)==1)lower.limits=rep(lower.limits,nvars)else stop("Require length 1 or nvars lower.limits")
  }
  else lower.limits=lower.limits[seq(nvars)]
  if(length(upper.limits)<nvars){
    if(length(upper.limits)==1)upper.limits=rep(upper.limits,nvars)else stop("Require length 1 or nvars upper.limits")
  }
  else upper.limits=upper.limits[seq(nvars)]
  cl=rbind(lower.limits,upper.limits)
  if(any(cl==0)){
    ###Bounds of zero can mess with the lambda sequence and fdev; ie nothing happens and if fdev is not
    ###zero, the path can stop
    fdev=glmnet.control()$fdev
    if(fdev!=0) {
      glmnet.control(fdev=0)
      on.exit(glmnet.control(fdev=fdev))
    }
  }
  storage.mode(cl)="double"
  ### end check on limits
  
  isd=as.integer(standardize)
  intr=as.integer(intercept)
  if(!missing(intercept)&&family=="cox")warning("Cox model has no intercept")
  jsd=as.integer(standardize.response)
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
  if(inherits(x,"sparseMatrix")){##Sparse case
    is.sparse=TRUE
    x=as(x,"CsparseMatrix")
    x=as(x,"dgCMatrix")
    ix=as.integer(x@p+1)
    jx=as.integer(x@i+1)
    x=as.double(x@x)
  }


  kopt=switch(match.arg(type.logistic),
   "Newton"=0,#This means to use the exact Hessian
    "modified.Newton"=1 # Use the upper bound
    )
  if(family=="multinomial"){
      type.multinomial=match.arg(type.multinomial)
      if(type.multinomial=="grouped")kopt=2 #overrules previous kopt
    }
  kopt=as.integer(kopt)

  fit=switch(family,
    "gaussian"=elnet(x,is.sparse,ix,jx,y,weights,offset,type.gaussian,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,vnames,maxit),
    "poisson"=fishnet(x,is.sparse,ix,jx,y,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,vnames,maxit),
    "binomial"=lognet(x,is.sparse,ix,jx,y,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,vnames,maxit,kopt,family),
    "multinomial"=lognet(x,is.sparse,ix,jx,y,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,vnames,maxit,kopt,family),
    "cox"=coxnet(x,is.sparse,ix,jx,y,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,vnames,maxit),
    "mgaussian"=mrelnet(x,is.sparse,ix,jx,y,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,jsd,intr,vnames,maxit)
    )
  if(is.null(lambda))fit$lambda=fix.lam(fit$lambda)##first lambda is infinity; changed to entry point
fit$call=this.call
  fit$nobs=nobs
  class(fit)=c(class(fit),"glmnet")
  fit
}

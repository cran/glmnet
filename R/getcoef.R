getcoef=function(fit,nvars,nx,vnames){
  lmu=fit$lmu
  if(lmu<1){
    ## changed this to a warning message, and return an empty model
    warning("an empty model has been returned; probably a convergence issue")
    coefob=list(a0=fit$a0,beta=zeromat(nvars,as.integer(1),vnames,"s0"),df=0,dim=c(nvars,1),lambda=Inf)
    return(coefob)
  }
  nin=fit$nin[seq(lmu)]
  ninmax=max(nin)
  lam=fit$alm[seq(lmu)]

  stepnames=paste("s",seq(lmu)-1,sep="")
  dd=c(nvars,lmu)
  if(ninmax>0){
           ca=matrix(fit$ca[seq(nx*lmu)],nx,lmu)[seq(ninmax),,drop=FALSE]
           df=apply(abs(ca)>0,2,sum)
           ja=fit$ia[seq(ninmax)]
####confusing but too hard to change
###glmnet builds a list of ever active variables which is nondecreasing
###Since ca was initialized to zero, no harm is done in passing a square matrix
###to new(); then when we do a drop0 it makes it really sparse
           oja=order(ja)
           ja=rep(ja[oja],lmu)
           ia=cumsum(c(1,rep(ninmax,lmu)))
           beta=drop0(new("dgCMatrix",Dim=dd,Dimnames=list(vnames,stepnames),x=as.vector(ca[oja,]),p=as.integer(ia-1),i=as.integer(ja-1)))
         }else {
           beta = zeromat(nvars,lmu,vnames,stepnames)
           df=rep(0,lmu)
         }
  a0=fit$a0
  if(!is.null(a0)){#for Cox model
    a0=a0[seq(lmu)]
  names(a0)=stepnames
  }
  list(a0=a0,beta=beta,df=df,dim=dd,lambda=lam)
}

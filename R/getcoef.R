getcoef=function(fit,nvars,nx,vnames){
  lmu=fit$lmu
  nin=fit$nin[seq(lmu)]
  ninmax=max(nin)
  lam=fit$alm[seq(lmu)]

  stepnames=paste("s",seq(lmu)-1,sep="")
  dd=c(nvars,lmu)
  if(ninmax>0){
           ca=matrix(fit$ca[seq(nx*lmu)],nx,lmu)[seq(ninmax),,drop=FALSE]
           df=apply(abs(ca)>0,2,sum)
           ja=fit$ia[seq(ninmax)]#confusing but too hard to change
           oja=order(ja)
           ja=rep(ja[oja],lmu)
           ia=cumsum(c(1,rep(ninmax,lmu)))
           beta=new("dgCMatrix",Dim=dd,Dimnames=list(vnames,stepnames),x=as.vector(ca[oja,]),p=as.integer(ia-1),i=as.integer(ja-1))
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

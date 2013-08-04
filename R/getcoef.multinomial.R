getcoef.multinomial=function(fit,nvars,nx,vnames,nc,classnames,center.intercept=TRUE){
  lmu=fit$lmu
  nin=fit$nin[seq(lmu)]
  ninmax=max(nin)
  lam=fit$alm[seq(lmu)]
  stepnames=paste("s",seq(lmu)-1,sep="")
  beta.list=as.list(seq(nc))
  names(beta.list)=classnames
    a0=matrix(fit$a0[seq(lmu*nc)],nc,lmu,dimnames=list(classnames,stepnames))
      if(center.intercept){
        a0=scale(a0,TRUE,FALSE)
        attr(a0,"scaled:center")=NULL
      }
      dfmat=a0
      dd=c(nvars,lmu)
      if(ninmax>0){
        ca=array(fit$ca[seq(nx*lmu*nc)],c(nx,nc,lmu))[seq(ninmax),,,drop=FALSE]
        ja=fit$ia[seq(ninmax)]#confusing but too hard to change
        oja=order(ja)
        ja=rep(ja[oja],lmu)
        ia=cumsum(c(1,rep(ninmax,lmu)))
        df=apply(abs(ca)>0,c(1,3),any)
        df=apply(df,2,sum)
        for(k in seq(nc)){
          cak=ca[oja,k, ,drop=FALSE]
          dfmat[k,]=apply(abs(cak)>0,3,sum)
          beta=new("dgCMatrix",Dim=dd,Dimnames=list(vnames,stepnames),x=as.vector(cak),p=as.integer(ia-1),i=as.integer(ja-1))
          beta.list[[k]]=drop0(beta)
        }
   } else{
    for (k in seq(nc)) {
      dfmat[k, ] = rep(0,lmu)
      beta.list[[k]] = zeromat(nvars,lmu,vnames,stepnames)
    }
  }
       
  list(a0=a0,beta=beta.list,dfmat=dfmat,df=df,dim=dd,lambda=lam)
}

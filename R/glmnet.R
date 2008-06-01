glmnet=function(x,y,family=c("gaussian","binomial","multinomial"),weights,alpha=1.0,nlambda=100,lambda.min=ifelse(nobs<nvars,5e-2,1e-4),lambda,standardize=TRUE,thresh=1e-4,dfmax=nvars+1,pmax=min(dfmax*1.2,nvars),exclude,penalty.factor=rep(1,nvars),maxit=100,HessianExact=FALSE,type=c("covariance","naive")){
  family=match.arg(family)
  this.call=match.call()
  nlam=as.integer(nlambda)
  np=dim(x)
  nobs=as.integer(np[1])
  if(missing(weights))weights=rep(1,nobs)
  nvars=as.integer(np[2])
  vnames=colnames(x)
  if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")
  if(family %in% c("binomial","multinomial")){
    nc=dim(y)
    maxit=as.integer(maxit)
    kopt=as.integer(HessianExact)
    if(is.null(nc)){
      ## Need to construct a y matrix, and include the weights
      y=as.factor(y)
      ntab=table(y)
      classnames=names(ntab)
      nc=as.integer(length(ntab))
      y=diag(nc)[as.numeric(y),]
    }
    else{
      noo=nc[1]
      if(noo!=nobs)stop("x and y have different number of rows")
      nc=as.integer(nc[2])
      classnames=colnames(y)
    }
    if(family=="binomial"){
      if(nc>2)stop("More than two classes; use multinomial family instead")
      nc=as.integer(1) # for calling multnet
    }
    if(!missing(weights))y=y*weights
       ### Compute the null deviance
    prior=apply(y,2,sum)
    sumw=sum(y)
    prior=prior/sumw
    nulldev= -2*sum(y*outer(rep(1,nobs),log(prior),"*"))/sumw

  }
else
     {
       weights=as.double(weights)
       ### compute the null deviance
       ybar=weighted.mean(y,weights)
       nulldev=weighted.mean( (y-ybar)^2,weights)
       type=match.arg(type)

       ka=as.integer(switch(type,
         covariance=1,
         naive=2,
         ))
     }
 storage.mode(y)="double"

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
    if(missing(lambda)){
      if(lambda.min>=1)stop("lambda.min should be less than 1")
      flmin=as.double(lambda.min)
      ulam=double(1)
    }
    else{
      flmin=as.double(1)    
      if(any(lambda<0))stop("lambdas should be non-negative")
      ulam=as.double(rev(sort(lambda)))
      nlam=as.integer(length(lambda))
    }
    is.sparse=FALSE
    if(class(x)=="dgCMatrix"){##Sparse case
      is.sparse=TRUE
      ix=as.integer(x@p+1)
      jx=as.integer(x@i+1)
      x=as.double(x@x)

      if(family=="gaussian")
              fit=.Fortran("spelnet",
        ka,parm=alpha,nobs,nvars,x,ix,jx,y,weights,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,
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

        else
                fit=.Fortran("splognet",
        parm=alpha,nobs,nvars,nc,x,ix,jx,y,jd,vp,ne=ne,nx,nlam,flmin,ulam,thresh,isd,maxit,kopt,
        lmu=integer(1),
        a0=double(nlam*nc),
        ca=double(nx*nlam*nc),
        ia=integer(nx),
        nin=integer(nlam),
        dev=double(nlam),
        alm=double(nlam),
        nlp=integer(1),
        jerr=integer(1),PACKAGE="glmnet"
        )

    }## end of sparse case
  else
     {
       ##regular nonsparse case
       if(family=="gaussian")
         fit=.Fortran("elnet",
          ka,parm=alpha,nobs,nvars,as.double(x),y,weights,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,
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
     else
               fit=.Fortran("lognet",
          parm=alpha,nobs,nvars,nc,as.double(x),y,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,maxit,kopt,
          lmu=integer(1),
          a0=double(nlam*nc),
          ca=double(nx*nlam*nc),
          ia=integer(nx),
          nin=integer(nlam),
          dev=double(nlam),
          alm=double(nlam),
          nlp=integer(1),
          jerr=integer(1),PACKAGE="glmnet"
          )
     }## end regular nonsparse case
     lmu=fit$lmu
     nin=fit$nin[seq(lmu)]
     ninmax=max(nin)
     lam=fit$alm[seq(lmu)]
     if(missing(lambda))lam=fix.lam(lam)##first lambda is infinity; changed to entry point
  stepnames=paste("s",seq(lmu)-1,sep="")

     errmsg=jerr(fit$jerr,maxit,pmax)### error messages from fortran
     switch(paste(errmsg$n),
            "1"=stop(errmsg$msg,call.=FALSE),
            "-1"=warning(errmsg$msg,call.=FALSE)
            )
if(family=="multinomial"){
      beta.list=as.list(seq(nc))
      names(beta.list)=classnames
      a0=matrix(fit$a0[seq(lmu*nc)],nc,lmu,dimnames=list(classnames,stepnames))
      a0=scale(a0,TRUE,FALSE)
      attr(a0,"scaled:center")=NULL
    dfmat=a0
    ca=array(fit$ca[seq(nx*lmu*nc)],c(nx,nc,lmu))[seq(ninmax),,,drop=FALSE]
    ia=fit$ia[seq(ninmax)]
    ja=rep(fit$ia[seq(ninmax)],lmu)
    ia=cumsum(c(1,rep(ninmax,lmu)))
    df=apply(abs(ca)>0,c(1,3),any)
    df=apply(df,2,sum)
    dd=c(nvars,lmu)
     for(k in seq(nc)){
      cak=ca[,k, ]
      dfmat[k,]=apply(abs(cak)>0,2,sum)
      beta=new("dgCMatrix",Dim=dd,Dimnames=list(vnames,stepnames),x=as.vector(cak),p=as.integer(ia-1),i=as.integer(ja-1))
      beta.list[[k]]=beta
    }
    outlist=list(a0=a0,beta=beta.list,
         dev=fit$dev[seq(lmu)],nulldev=nulldev,dfmat=dfmat,df=df,
         lambda=lam,npasses=fit$nlp,jerr=fit$jerr,dim=dd,call=this.call)
    class(outlist)=c("glmnet","multnet")
    }
     else{
           ca=matrix(fit$ca[seq(nx*lmu)],nx,lmu)[seq(ninmax),]
           df=apply(abs(ca)>0,2,sum)
           ia=fit$ia[seq(ninmax)]
           ja=rep(fit$ia[seq(ninmax)],lmu)
           ia=cumsum(c(1,rep(ninmax,lmu)))
           dd=c(nvars,lmu)
           beta=new("dgCMatrix",Dim=dd,Dimnames=list(vnames,stepnames),x=as.vector(ca),p=as.integer(ia-1),i=as.integer(ja-1))

      if(family=="binomial"){
          a0=-fit$a0[seq(lmu)]
          names(a0)=stepnames
           outlist=list(a0=a0,beta=-beta,#sign flips make 2 target class
         dev=fit$dev[seq(lmu)],nulldev=nulldev,df=df,
         lambda=lam,npasses=fit$nlp,jerr=fit$jerr,dim=dd,call=this.call)
    class(outlist)=c("glmnet","lognet")
         }
           else
             {
               a0=fit$a0[seq(lmu)]
               names(a0)=stepnames
                  outlist= list(a0=a0,beta=beta,dev=fit$rsq[seq(lmu)],nulldev=nulldev,df=df,
         lambda=lam,npasses=fit$nlp,jerr=fit$jerr,dim=dd,call=this.call)
    class(outlist)=c("glmnet","elnet")
                }
         }
     outlist
   }

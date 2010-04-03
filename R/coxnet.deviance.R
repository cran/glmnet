coxnet.deviance=function(pred=NULL,y,x=0,offset=pred,weights=NULL,beta=NULL){
  storage.mode(x)="double"
  y=response.coxnet(y)
  ty=y$time
  tevent=y$event
  nobs=as.integer(length(ty))
  nvars=as.integer(ncol(x))
  if(is.null(weights))weights=rep(1.0,nobs)
  else weights=as.double(weights)
  if(is.null(offset))offset=rep(0.0,nobs)
  else offset=as.double(offset)
  if(is.null(beta)){
    beta=double(0)
    nvec=as.integer(1)
    nvars=as.integer(0)
  }
  else{
    beta=as.matrix(beta)
    nvec=as.integer(ncol(beta))
  }
  fit=.Fortran("loglike",
    nobs,nvars,x,ty,tevent,offset,weights,nvec,beta,
    flog=double(nvec),
    jerr=integer(1),
    PACKAGE="glmnet")
  if(fit$jerr!=0){
  errmsg=jerr(fit$jerr,maxit=0,pmax=0,family="cox")
  if(errmsg$fatal)stop(errmsg$msg,call.=FALSE)
  else warning(errmsg$msg,call.=FALSE)
}
  
-2*fit$flog
}

  
  

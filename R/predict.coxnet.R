predict.coxnet=function(object,newx,s=NULL,type=c("link","response","coefficients","nonzero"),exact=FALSE,offset,...){
  type=match.arg(type)
  ###coxnet has no intercept, so we treat it separately
  if(missing(newx)){
    if(!match(type,c("coefficients","nonzero"),FALSE))stop("You need to supply a value for 'newx'")
  }
  if(exact&&(!is.null(s))){
    lambda=object$lambda
    which=match(s,lambda,FALSE)
    if(!all(which>0)){
      lambda=unique(rev(sort(c(s,lambda))))
      object=update(object,lambda=lambda)
    }
  }

  nbeta=object$beta
   if(!is.null(s)){
    vnames=dimnames(nbeta)[[1]]
    dimnames(nbeta)=list(NULL,NULL)
    lambda=object$lambda
    lamlist=lambda.interp(lambda,s)
    nbeta=nbeta[,lamlist$left,drop=FALSE]*lamlist$frac +nbeta[,lamlist$right,drop=FALSE]*(1-lamlist$frac)
    dimnames(nbeta)=list(vnames,paste(seq(along=s)))
  }
  if(type=="coefficients")return(nbeta)
  if(type=="nonzero")return(nonzeroCoef(nbeta,bystep=TRUE))
  nfit=as.matrix(newx%*%nbeta)
  if(object$offset){
    if(missing(offset))stop("No offset provided for prediction, yet used in fit of glmnet",call.=FALSE)
    nfit=nfit+array(offset,dim=dim(nfit))
  }
  switch(type,
         response=exp(nfit),
         link=nfit
         )
}  

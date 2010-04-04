predict.fishnet=function(object,newx,s=NULL,type=c("link","response","coefficients","nonzero"),exact=FALSE,offset,...){
  type=match.arg(type)
  a0=t(as.matrix(object$a0))
  rownames(a0)="(Intercept)"
  nbeta=rbind2(a0,object$beta)
  if(!is.null(s)){
    vnames=dimnames(nbeta)[[1]]
    dimnames(nbeta)=list(NULL,NULL)
    lambda=object$lambda
    lamlist=lambda.interp(lambda,s)
    nbeta=nbeta[,lamlist$left,drop=FALSE]*lamlist$frac +nbeta[,lamlist$right,drop=FALSE]*(1-lamlist$frac)
    dimnames(nbeta)=list(vnames,paste(seq(along=s)))
  }
  if(type=="coefficients")return(nbeta)
  if(type=="nonzero")return(nonzeroCoef(nbeta[-1,,drop=FALSE],bystep=TRUE))
  nfit=as.matrix(cbind2(1,newx)%*%nbeta)
  if(object$offset){
    if(missing(offset))stop("No offset provided for prediction, yet used in fit of glmnet",call.=FALSE)
    nfit=nfit+array(offset,dim=dim(nfit))
  }
  switch(type,
         response=exp(nfit),
         link=nfit
         )
}  

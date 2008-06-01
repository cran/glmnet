predict.lognet=function(object,newx,s=object$lambda,type=c("link","response","coefficients","class","nonzero"),exact=FALSE,...){
  type=match.arg(type)
  a0=t(as.matrix(object$a0))
  rownames(a0)="(Intercept)"
  nbeta=rbind2(object$a0,object$beta)
  if(!missing(s)){
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
  switch(type,
         response={
           pp=exp(-nfit)
           1/(1+pp)
         },
         link=nfit,
         class=ifelse(nfit>0,1,2)
         )
}  
    

predict.elnet=function(object,newx,s=object$lambda,type=c("link","response","coefficients","nonzero"),exact=FALSE,...){
  type=match.arg(type)
  a0=t(as.matrix(object$a0))
  rownames(a0)="(Intercept)"
  nbeta=rbind2(a0,object$beta)
  if(!missing(s)){
    vnames=dimnames(nbeta)[[1]]
    dimnames(nbeta)=list(NULL,NULL)
    lambda=object$lambda
    lamlist=lambda.interp(lambda,s)
    nbeta=nbeta[,lamlist$left,drop=FALSE]*lamlist$frac +nbeta[,lamlist$right,drop=FALSE]*(1-lamlist$frac)
    dimnames(nbeta)=list(vnames,paste(seq(along=s)))
  }
  switch(type,
         "coefficients"=nbeta,
         "link"=as.matrix(cbind2(1,newx)%*%nbeta),
         "response"=as.matrix(cbind2(1,newx)%*%nbeta),
         "nonzero"=nonzeroCoef(nbeta[-1,,drop=FALSE],bystep=TRUE)
         )
}  
    

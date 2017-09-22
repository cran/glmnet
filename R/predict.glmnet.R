predict.glmnet=function(object,newx,s=NULL,type=c("link","response","coefficients","nonzero","class"),exact=FALSE,newoffset,...){
 type=match.arg(type)
  if(missing(newx)){
    if(!match(type,c("coefficients","nonzero"),FALSE))stop("You need to supply a value for 'newx'")
     }
  if(exact&&(!is.null(s))){
###we augment the lambda sequence with the new values, if they are different,and refit the model using update
    lambda=object$lambda
    which=match(s,lambda,FALSE)
    if(!all(which>0)){
        lambda=unique(rev(sort(c(s,lambda))))
        check_dots(object,...)# This fails if you don't supply the crucial arguments
        object=update(object,lambda=lambda,...)
    }
  }
  a0=t(as.matrix(object$a0))
  rownames(a0)="(Intercept)"
  nbeta=methods::rbind2(a0,object$beta)#was rbind2
  if(!is.null(s)){
    vnames=dimnames(nbeta)[[1]]
    dimnames(nbeta)=list(NULL,NULL)
    lambda=object$lambda
    lamlist=lambda.interp(lambda,s)

    nbeta=nbeta[,lamlist$left,drop=FALSE]%*%Diagonal(x=lamlist$frac) +nbeta[,lamlist$right,drop=FALSE]%*%Diagonal(x=1-lamlist$frac)
    dimnames(nbeta)=list(vnames,paste(seq(along=s)))
  }
  if(type=="coefficients")return(nbeta)
  if(type=="nonzero")return(nonzeroCoef(nbeta[-1,,drop=FALSE],bystep=TRUE))
  ###Check on newx
 if(inherits(newx, "sparseMatrix"))newx=as(newx,"dgCMatrix")
  nfit=as.matrix(cbind2(1,newx)%*%nbeta)
   if(object$offset){
    if(missing(newoffset))stop("No newoffset provided for prediction, yet offset used in fit of glmnet",call.=FALSE)
    if(is.matrix(newoffset)&&inherits(object,"lognet")&&dim(newoffset)[[2]]==2)newoffset=newoffset[,2]
    nfit=nfit+array(newoffset,dim=dim(nfit))
  }
nfit
  }

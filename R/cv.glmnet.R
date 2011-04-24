cv.glmnet=function(x,y,weights,offset=NULL,lambda=NULL,type.measure=c("mse","deviance","class","auc","mae"),...,nfolds=10,foldid,grouped=TRUE){
  if(missing(type.measure))type.measure="default"
  else type.measure=match.arg(type.measure)
  N=nrow(x)
  if(missing(weights))weights=rep(1.0,N)else weights=as.double(weights)
###Fit the model once to get dimensions etc of output
  y=drop(y) # we dont like matrix responses unless we need them
  glmnet.object=glmnet(x,y,weights=weights,offset=offset,lambda=lambda,...)
  is.offset=glmnet.object$offset
  lambda=glmnet.object$lambda
  if(inherits(glmnet.object,"multnet")){
    nz=predict(glmnet.object,type="nonzero")
    nz=sapply(nz,function(x)sapply(x,length))
    nz=ceiling(apply(nz,1,median))
  }
  else  nz=sapply(predict(glmnet.object,type="nonzero"),length)
  if(missing(foldid)) foldid=sample(rep(seq(nfolds),length=N)) else nfolds=max(foldid)
  if(nfolds<3)stop("nfolds must be bigger than 3; nfolds=10 recommended")
   outlist=as.list(seq(nfolds))
###Now fit the nfold models and store them
  for(i in seq(nfolds)){
    which=foldid==i
    if(is.matrix(y))y_sub=y[!which,]else y_sub=y[!which]
    if(is.offset)offset_sub=as.matrix(offset)[!which,]
    else offset_sub=NULL
    outlist[[i]]=glmnet(x[!which,,drop=FALSE],y_sub,lambda=lambda,offset=offset_sub,weights=weights[!which],...)
  }
  ###What to do depends on the type.measure and the model fit
  fun=paste("cv",class(glmnet.object)[[1]],sep=".")
  cvstuff=do.call(fun,list(outlist,lambda,x,y,weights,offset,foldid,type.measure,grouped))
  cvm=cvstuff$cvm
  cvsd=cvstuff$cvsd
  cvname=cvstuff$name
  
out=list(lambda=lambda,cvm=cvm,cvsd=cvsd,cvup=cvm+cvsd,cvlo=cvm-cvsd,nzero=nz,name=cvname,glmnet.fit=glmnet.object)
  lamin=if(type.measure=="auc")getmin(lambda,-cvm,cvsd)
  else getmin(lambda,cvm,cvsd)
  obj=c(out,as.list(lamin))
  class(obj)="cv.glmnet"
  obj
}

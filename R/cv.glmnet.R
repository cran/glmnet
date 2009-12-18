
cv.glmnet=function(x,y,...,nfolds=10,foldid,type=c("response","class")){
  type=match.arg(type)
###Fit the model once to get dimensions etc of output
  glmnet.object=glmnet(x,y,...)
  lambda=glmnet.object$lambda
  N=nrow(x)
  if(missing(foldid)) foldid=sample(rep(seq(nfolds),length=N)) else nfolds=max(foldid)
  predmat=predict(glmnet.object,x,lambda=lambda)
  nz=sapply(predict(glmnet.object,type="nonzero"),length)

  for(i in seq(nfolds)){
    which=foldid==i
    cvfit=glmnet(x[!which,],y[!which],lambda=lambda,...)
    predmat[which,]=predict(cvfit,x[which,],type=type)
  }
  ###What to do depends on the type and the model fit
 cvstuff= switch(class(glmnet.object)[[2]],
         elnet=cvelnet(y,predmat),
         lognet=cvlognet(y,predmat,type),
         multnet=stop("not yet implemented")
         )
 cvm=cvstuff$cvm
 cvsd=cvstuff$cvsd
  cvname=cvstuff$name
  
out=list(lambda=lambda,cvm=cvm,cvsd=cvsd,cvup=cvm+cvsd,cvlo=cvm-cvsd,nzero=nz,name=cvname)
  lamin=getmin(lambda,cvm,cvsd)
  obj=c(out,as.list(lamin))
  class(obj)="cv.glmnet"
  obj
}

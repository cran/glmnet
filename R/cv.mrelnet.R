cv.mrelnet=function(outlist,lambda,x,y,weights,offset,foldid,type.measure,grouped,keep=FALSE){
  typenames=c(deviance="Mean-Squared Error",mse="Mean-Squared Error",mae="Mean Absolute Error")
  if(type.measure=="default")type.measure="mse"
  if(!match(type.measure,c("mse","mae","deviance"),FALSE)){
    warning("Only 'mse', 'deviance' or 'mae'  available for multiresponse Gaussian models; 'mse' used")
    type.measure="mse"
  }
  typename=typenames[type.measure]
  if(type.measure=="deviance")type.measure="mse"
  ndim=dim(y)
  nc=ndim[2]
  nobs=ndim[1]
  if(!is.null(offset))y=y-drop(offset)
  predmat=array(NA,c(nobs,nc,length(lambda)))
    nfolds=max(foldid)
    nlams=double(nfolds)
    for(i in seq(nfolds)){
      which=foldid==i
      fitobj=outlist[[i]]
      fitobj$offset=FALSE
      preds=predict(fitobj,x[which,,drop=FALSE])
      nlami=length(outlist[[i]]$lambda)
       predmat[which,,seq(nlami)]=preds
      nlams[i]=nlami
    }

  N=nobs - apply(is.na(predmat[,1,]),2,sum)
  bigY=array(y,dim(predmat))
  cvraw=switch(type.measure,
    "mse"=apply((bigY-predmat)^2,c(1,3),sum),
    "mae"=apply(abs(bigY-predmat),c(1,3),sum)
    )
   if( (nobs/nfolds <3)&&grouped){
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",call.=FALSE)
    grouped=FALSE
  }
 if(grouped){
   cvob=cvcompute(cvraw,weights,foldid,nlams)
  cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
 }

  cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  out=list(cvm=cvm,cvsd=cvsd,name=typename)
  if(keep)out$fit.preval=predmat
  out

}

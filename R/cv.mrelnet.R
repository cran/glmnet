cv.mrelnet=function(outlist,lambda,x,y,weights,offset,foldid,type.measure,grouped,keep=FALSE){
  ndim=dim(y)
  nc=ndim[2]
  nobs=ndim[1]
  if(!is.null(offset))y=y-drop(offset)
  ##We dont want to extrapolate lambdas on the small side
  mlami=max(sapply(outlist,function(obj)min(obj$lambda)))
  which_lam=lambda >= mlami

  predmat=array(NA,c(nobs,nc,length(lambda)))
    nfolds=max(foldid)
    nlams=double(nfolds)
    for(i in seq(nfolds)){
      which=foldid==i
      fitobj=outlist[[i]]
      fitobj$offset=FALSE
      preds=predict(fitobj,x[which,,drop=FALSE], s=lambda[which_lam])
      nlami=sum(which_lam)
       predmat[which,,seq(nlami)]=preds
      nlams[i]=nlami
    }

  N=nobs - apply(is.na(predmat[,1,]),2,sum)
  bigY=array(y,dim(predmat))
  cvraw=switch(type.measure,
    "mse"=apply((bigY-predmat)^2,c(1,3),sum),
    "deviance"=apply((bigY-predmat)^2,c(1,3),sum),
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
  out=list(cvm=cvm,cvsd=cvsd,type.measure=type.measure)
  if(keep)out$fit.preval=predmat
  out

}

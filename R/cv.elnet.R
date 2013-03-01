cv.elnet=function(outlist,lambda,x,y,weights,offset,foldid,type.measure,grouped,keep=FALSE){
  typenames=c(deviance="Mean-Squared Error",mse="Mean-Squared Error",mae="Mean Absolute Error")
  if(type.measure=="default")type.measure="mse"
  if(!match(type.measure,c("mse","mae","deviance"),FALSE)){
    warning("Only 'mse', 'deviance' or 'mae'  available for Gaussian models; 'mse' used")
    type.measure="mse"
  }
     if(!is.null(offset))y=y-drop(offset)
     predmat=matrix(NA,length(y),length(lambda))
    nfolds=max(foldid)
    nlams=double(nfolds)
    for(i in seq(nfolds)){
      which=foldid==i
      fitobj=outlist[[i]]
      fitobj$offset=FALSE
      preds=predict(fitobj,x[which,,drop=FALSE])
      nlami=length(outlist[[i]]$lambda)
       predmat[which,seq(nlami)]=preds
      nlams[i]=nlami
    }

  N=length(y) - apply(is.na(predmat),2,sum)
  cvraw=switch(type.measure,
    "mse"=(y-predmat)^2,
    "deviance"=(y-predmat)^2,
    "mae"=abs(y-predmat)
    )
   if( (length(y)/nfolds <3)&&grouped){
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",call.=FALSE)
    grouped=FALSE
  }
 if(grouped){
   cvob=cvcompute(cvraw,weights,foldid,nlams)
  cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
 }

  cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  out=list(cvm=cvm,cvsd=cvsd,name=typenames[type.measure])
  if(keep)out$fit.preval=predmat
  out
}

cv.fishnet=function(outlist,lambda,x,y,weights,offset,foldid,type.measure,grouped,keep=FALSE){
  typenames=c(mse="Mean-Squared Error",mae="Mean Absolute Error","deviance"="Poisson Deviance")
  if(type.measure=="default")type.measure="deviance"
  if(!match(type.measure,c("mse","mae","deviance"),FALSE)){
    warning("Only 'deviance', 'mse' or 'mae'  available for Poisson models; 'deviance' used")
    type.measure="deviance"
  }
     if(!is.null(offset)){
       is.offset=TRUE
       offset=drop(offset)
     }else is.offset=FALSE

  devi=function(y,eta){
    deveta=y*eta-exp(eta)
    devy=y*log(y)-y
    devy[y==0]=0
    2*(devy-deveta)
  }
     predmat=matrix(NA,length(y),length(lambda))
    nfolds=max(foldid)
  nlams=double(nfolds)
    for(i in seq(nfolds)){
      which=foldid==i
      fitobj=outlist[[i]]
      if(is.offset)off_sub=offset[which]
      preds=predict(fitobj,x[which,,drop=FALSE],offset=off_sub)
      nlami=length(outlist[[i]]$lambda)
      predmat[which,seq(nlami)]=preds
      nlams[i]=nlami
    }

  N=length(y) - apply(is.na(predmat),2,sum)
  cvraw=switch(type.measure,
    "mse"=(y-exp(predmat))^2,
    "mae"=abs(y-exp(predmat)),
    "deviance"=devi(y,predmat)
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

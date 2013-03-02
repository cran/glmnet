cv.lognet=function(outlist,lambda,x,y,weights,offset,foldid,type.measure,grouped,keep=FALSE){
  typenames=c(mse="Mean-Squared Error",mae="Mean Absolute Error",deviance="Binomial Deviance",auc="AUC",class="Misclassification Error")
  if(type.measure=="default")type.measure="deviance"
  if(!match(type.measure,c("mse","mae","deviance","auc","class"),FALSE)){
    warning("Only 'deviance', 'class', 'auc', 'mse' or 'mae'  available for binomial models; 'deviance' used")
    type.measure="deviance"
  }

###These are hard coded in the Fortran, so we do that here too
  prob_min=1e-5
  prob_max=1-prob_min
  ###Turn y into a matrix
  nc = dim(y)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
    }
  N=nrow(y)
  nfolds=max(foldid)
  if( (N/nfolds <10)&&type.measure=="auc"){
    warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
    type.measure="deviance"
  }
  if( (N/nfolds <3)&&grouped){
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",call.=FALSE)
    grouped=FALSE
  }


  if(!is.null(offset)){
    is.offset=TRUE
    offset=drop(offset)
  }else is.offset=FALSE
  predmat=matrix(NA,nrow(y),length(lambda))
   nlams=double(nfolds)
  for(i in seq(nfolds)){
      which=foldid==i
      fitobj=outlist[[i]]
      if(is.offset)off_sub=offset[which]
      preds=predict(fitobj,x[which,,drop=FALSE],offset=off_sub,type="response")
      nlami=length(outlist[[i]]$lambda)
      predmat[which,seq(nlami)]=preds
      nlams[i]=nlami
    }
   ###If auc we behave differently
  if(type.measure=="auc"){
    cvraw=matrix(NA,nfolds,length(lambda))
    good=matrix(0,nfolds,length(lambda))
    for(i in seq(nfolds)){
      good[i,seq(nlams[i])]=1
      which=foldid==i
      for(j in seq(nlams[i])){
        cvraw[i,j]=auc.mat(y[which,],predmat[which,j],weights[which])
      }
    }
    N=apply(good,2,sum)
    weights=tapply(weights,foldid,sum)
  }
   else{
    ##extract weights and normalize to sum to 1
    ywt=apply(y,1,sum)
    y=y/ywt
    weights=weights*ywt

    N=nrow(y) - apply(is.na(predmat),2,sum)
    cvraw=switch(type.measure,
    "mse"=(y[,1]-(1-predmat))^2 +(y[,2]-predmat)^2,
    "mae"=abs(y[,1]-(1-predmat)) +abs(y[,2]-predmat),
    "deviance"= {
      predmat=pmin(pmax(predmat,prob_min),prob_max)
      lp=y[,1]*log(1-predmat)+y[,2]*log(predmat)
      ly=log(y)
      ly[y==0]=0
      ly=drop((y*ly)%*%c(1,1))
      2*(ly-lp)
    },
    "class"=y[,1]*(predmat>.5) +y[,2]*(predmat<=.5)
    )
 if(grouped){
   cvob=cvcompute(cvraw,weights,foldid,nlams)
  cvraw=cvob$cvraw;weights=cvob$weights;N=cvob$N
 }
  }
   cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  out=list(cvm=cvm,cvsd=cvsd,name=typenames[type.measure])
  if(keep)out$fit.preval=predmat
  out

}

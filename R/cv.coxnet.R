cv.coxnet=function(outlist,lambda,x,y,weights,offset,foldid,type,grouped){
  typenames=c("deviance"="Partial Likelihood Deviance")
  if(type=="default")type="deviance"
  if(!match(type,c("deviance"),FALSE)){
    warning("Only 'deviance'  available for Cox models; 'deviance' used")
    type="deviance"
  }
     if(!is.null(offset)){
       is.offset=TRUE
       offset=drop(offset)
     }else is.offset=FALSE


    nfolds=max(foldid)
  cvraw=matrix(NA,nfolds,length(lambda))
  for(i in seq(nfolds)){
    which=foldid==i
    fitobj=outlist[[i]]
    coefmat=predict(fitobj,type="coeff")
    if(grouped){
      plfull=coxnet.deviance(x=x,y=y,offset=offset,weights=weights,beta=coefmat)
      plminusk=coxnet.deviance(x=x[!which,],y=y[!which,],offset=offset[!which],weights=weights[!which],beta=coefmat)
      cvraw[i,seq(along=plfull)]=plfull-plminusk
    }
    else{
      plk=coxnet.deviance(x=x[which,],y=y[which,],offset=offset[which],weights=weights[which],beta=coefmat)
      cvraw[i,seq(along=plk)]=plk
    }
  }
  status=y[,"status"]
  N=nfolds - apply(is.na(cvraw),2,sum)
  weights=as.vector(tapply(weights*status,foldid,sum))
   cvraw=cvraw/weights

  cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  list(cvm=cvm,cvsd=cvsd,name=typenames[type])
}

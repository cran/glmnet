cvlognet=function(y,predmat,type=c("response","class")){
  N=dim(predmat)[[1]]
  ##assume for now that y is binary
  y=as.numeric(as.factor(y))
   cvraw=switch(type,
     response =-2*((y==2)*log(predmat)+(y==1)*log(1-predmat)),
     class=y!=predmat
     )
  cvm=apply(cvraw,2,mean)
  cvsd=sqrt(apply(cvraw,2,var)/N)
  cvname=switch(type,response="Deviance",class="Misclassification Error")
  list(cvm=cvm,cvsd=cvsd,name=cvname)
}

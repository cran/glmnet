nonzeroCoef=function(beta,bystep=FALSE){
  ##beta should be in "dgCMatrix" format
  beta=t(beta)
  which=diff(beta@p)
  which=seq(which)[which>0]
  if(bystep){
    nzel=function(x,which)if(any(x))which[x]else NULL
    beta=abs(as.matrix(beta[,which]))>0
    apply(beta,1,nzel,which)
  }
  else which
    
  }

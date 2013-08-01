jerr=function(n,maxit,pmax,family){
  if(n==0) list(n=0,fatal=FALSE,msg="")
  else {
    errlist=switch(family,
              "gaussian"=jerr.elnet(n,maxit,pmax),
               "binomial"=jerr.lognet(n,maxit,pmax),
              "multinomial"=jerr.lognet(n,maxit,pmax),
              "poisson"=jerr.fishnet(n,maxit,pmax),
              "cox"=jerr.coxnet(n,maxit,pmax),
              "mrelnet"=jerr.mrelnet(n,maxit,pmax)
      )
    names(errlist)=c("n","fatal","msg")
    errlist$msg=paste("from glmnet Fortran code (error code ",n, "); ",errlist$msg,sep="")
    errlist
  }
}
                  
                  
     

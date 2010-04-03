jerr=function(n,maxit,pmax,family){
  if(n==0) list(n=0,msg="",fatal=FALSE)
  else {
    errlist=switch(family,
              "gaussian"=jerr.elnet(n),
              "binomial"=jerr.lognet(n,maxit,pmax),
              "multinomial"=jerr.lognet(n,maxit,pmax),
              "poisson"=jerr.fishnet(n,maxit,pmax),
              "cox"=jerr.coxnet(n,maxit,pmax)
      )
    errlist$msg=paste("from glmnet Fortran code -", errlist$msg)
    errlist
  }
}
                  
                  
     

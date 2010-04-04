jerr.lognet=function(n,maxit,pmax){

  if(n>0){#fatal error
    outlist=jerr.elnet(n)
    if(outlist$msg!="Unknown error")return(outlist)
    if((8000<n) & (n<9000))msg=paste("Null probability for class",n-8000, "< 1.0e-5")
    else if((9000<n) & (n<10000))msg=paste("Null probability for class",n-9000, "> 1.0 - 1.0e-5")
    else msg="Unknown error"
    list(n,fatal=TRUE,msg)
  }
 if(n<0){ # non-fatal error
       if(n>-10000)msg=paste("Convergence for ",-n,"th lambda value not reached after maxit=",maxit," iterations; solutions for larger lambdas returned",sep="")
    if(n < -10000)msg=paste("Number of nonzero coefficients along the path exceeds pmax=",pmax, " at ",-n-10000,"th lambda value; solutions for larger lambdas returned",sep="")
list(n,fatal=FALSE,msg)
     }
}

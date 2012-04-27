jerr.elnet=function(n,maxit,pmax){
  if(n>0){#fatal error
    if(n<7777)msg="Memory allocation error; contact package maintainer"
    else if(n==7777)msg="All used predictors have zero variance"
    else if(n==10000)msg="All penalty factors are <= 0"
    else msg="Unknown error"
    list(n=n,fatal=TRUE,msg=msg)
  }
  else if(n<0){# non-fatal error
           if(n>-10000)msg=paste("Convergence for ",-n,"th lambda value not reached after maxit=",maxit," iterations; solutions for larger lambdas returned",sep="")
    if(n < -10000)msg=paste("Number of nonzero coefficients along the path exceeds pmax=",pmax, " at ",-n-10000,"th lambda value; solutions for larger lambdas returned",sep="")
list(n=n,fatal=FALSE,msg=msg)
     }

}

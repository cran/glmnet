jerr=function(n,maxit,pmax){
  if(n==0) msg=""
  if(n>0){#fatal error
    if(n<7777)msg="Memory allocation error; contact package maintainer"
    if(n==7777)msg="All used predictors have zero variance"
    if((8000<n) & (n<9000))msg=paste("Null probability for class",n-8000, "< 1.0e-5")
    if((9000<n) & (n<10000))msg=paste("Null probability for class",n-9000, "> 1.0 - 1.0e-5")
    if(n==10000)msg="All penalty factors are <= 0"
    n=1
  msg=paste("in glmnet fortran code -",msg)
  }
  if(n<0){#non fatal error
    if(n>-10000)msg=paste("Convergence for ",-n,"th lambda value not reached after maxit=",maxit," iterations; solutions for larger lambdas returned",sep="")
    if(n < -10000)msg=paste("Number of nonzero coefficients along the path exceeds pmax=",pmax, " at ",-n-10000,"th lambda value; solutions for larger lambdas returned",sep="")
    n=-1
  msg=paste("from glmnet fortran code -",msg)
  }
  list(n=n,msg=msg)
}
                  
                  
     

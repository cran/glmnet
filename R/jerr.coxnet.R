jerr.coxnet=function(n,maxit,pmax){

  if(n>0){#fatal error
    outlist=jerr.elnet(n)
    if(outlist$msg!="")return(outlist)
    if(n==8888)msg="All observations censored - cannot proceed"  
    else if(n==9999)msg="No positive observation weights"
else 	 msg=""
    list(n,fatal=TRUE,msg)
  }
 if(n<0){ # non-fatal error
   if(n==-1)msg=paste("Inner loop convergence not reached after",maxit,"iterations")
   else if(n==-2)msg=paste("Outer loop convergence not reached after",maxit,"iterations")
   else if(n==-3)msg=paste("Numerical error in cox model in glmnet")
   else if(n < -10000)msg=paste("Number of nonzero coefficients along the path exceeds pmax=",pmax, " at ",-n-10000,"th lambda value; solutions for larger lambdas returned",sep="")
   else msg="Unknown error"
    list(n,fatal=FALSE,msg)
 }
}

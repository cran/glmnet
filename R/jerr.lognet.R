jerr.lognet=function(n,maxit,pmax){
    outlist=jerr.elnet(n,maxit,pmax)
    if(n<  -20000)outlist$msg=paste("Max(p(1-p),1.0e-6 at ",-n-20000,"th value of lambda; solutions for larger values of lambda returned")
    if(outlist$msg!="Unknown error")return(outlist)
    if((8000<n) & (n<9000))msg=paste("Null probability for class",n-8000, "< 1.0e-5")
    else if((9000<n) & (n<10000))msg=paste("Null probability for class",n-9000, "> 1.0 - 1.0e-5")
    else msg="Unknown error"
    list(n=n,fatal=TRUE,msg=msg)
}

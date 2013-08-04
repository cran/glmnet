jerr.coxnet=function(n,maxit,pmax){

  if(n>0){#fatal error
    outlist=jerr.elnet(n)
    if(outlist$msg!="Unknown error")return(outlist)
    if(n==8888)msg="All observations censored - cannot proceed"  
    else if(n==9999)msg="No positive observation weights"
    else if(match(n,c(20000,30000),FALSE)) msg="Inititialization numerical error; probably too many censored observations"
    else 	 msg="Unknown error"
    list(n=n,fatal=TRUE,msg=msg)
  } else if(n<0){ # non-fatal error
    if(n<= -30000){
      msg=paste("Numerical error at ",-n-30000,"th lambda value; solutions for larger values of lambda returned",sep="")
      list(n=n,fatal=FALSE,msg=msg)
    }
    else
      jerr.elnet(n,maxit,pmax)
  }
}

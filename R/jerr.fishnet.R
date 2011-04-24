jerr.fishnet=function(n,maxit,pmax){
  outlist=jerr.elnet(n,maxit,pmax)
  if(outlist$msg!="Unknown error")return(outlist)
  if(n==8888)msg="Negative response values - should be counts"  
  else if(n==9999)msg="No positive observation weights"
  else 	 msg="Unknown error"
  list(n=n,fatal=TRUE,msg=msg)
}

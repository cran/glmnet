fix.lam=function(lam){
  llam=log(lam)
  lam[1]=exp(2*llam[2]-llam[3])
  lam
}

cvelnet=function(y,predmat){
  N=length(y)
   cvraw=(y-predmat)^2
  cvm=apply(cvraw,2,mean)
  cvsd=sqrt(apply(cvraw,2,var)/N)
  list(cvm=cvm,cvsd=cvsd,name="Mean Squared Error")
}

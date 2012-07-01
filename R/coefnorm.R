coefnorm=function(coeflist,q=1){
  ###coeflist comes from mrelnet or multnet
  nmat=coeflist[[1]]*0
  nclass=length(coeflist)
  for(i in seq(nclass)){
    nmat=nmat+abs(coeflist[[i]])^q
  }
  nmat^(1/q)
}

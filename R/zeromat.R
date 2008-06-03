zeromat=function(nvars,lmu,vnames,stepnames){
  ca=rep(0,lmu)
  ia=seq(lmu+1)
  ja=rep(1,lmu)
  dd=c(nvars,lmu)
  new("dgCMatrix", Dim = dd,
      Dimnames = list(vnames,stepnames),
      x = as.vector(ca),
      p = as.integer(ia - 1), i = as.integer(ja - 1))
}

cvstats=function(cvstuff,foldid,nfolds,lambda,nz,grouped,...){
    if (grouped){
      nlams=rep(dim(cvstuff$cvraw)[2],nfolds) ## All the same - this should go
      cvstuff= cvcompute(cvstuff, foldid, nlams)
}
  cvm=with(cvstuff,apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE))
  cvsd=with(cvstuff, sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                               w = weights, na.rm = TRUE)/(N - 1)))
   nas=is.na(cvsd)
  if(any(nas)){
    lambda=lambda[!nas]
    cvm=cvm[!nas]
    cvsd=cvsd[!nas]
    nz=nz[!nas]
  }
list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
                                                  cvsd, cvlo = cvm - cvsd, nzero = nz)
}

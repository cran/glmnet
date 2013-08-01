cvcompute=function(mat,weights,foldid,nlams){
  ###Computes the weighted mean and SD within folds, and hence the se of the mean
  wisum=tapply(weights,foldid,sum)
  nfolds=max(foldid)
  outmat=matrix(NA,nfolds,ncol(mat))
  good=matrix(0,nfolds,ncol(mat))
  mat[is.infinite(mat)]=NA#just in case some infinities crept in
  for(i in seq(nfolds)){
    mati=mat[foldid==i,,drop=FALSE]
    wi=weights[foldid==i]
    outmat[i,]=apply(mati,2,weighted.mean,w=wi,na.rm=TRUE)
    good[i,seq(nlams[i])]=1
    }
  N=apply(good,2,sum)
  list(cvraw=outmat,weights=wisum,N=N)
}

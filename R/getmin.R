getmin=function(lambda,cvm,cvsd){
  cvmin=min(cvm)
  idmin=cvm<=cvmin
  lambda.min=max(lambda[idmin])
  idmin=match(lambda.min,lambda)
  semin=(cvm+cvsd)[idmin]
  idmin=cvm<semin
  lambda.1se=max(lambda[idmin])
  list(lambda.min=lambda.min,lambda.1se=lambda.1se)
}
 

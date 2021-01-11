getOptcv.glmnet <-
    function (lambda, cvm, cvsd,cvname)
{
    if(match(cvname,c("AUC","C-index"),0))cvm=-cvm
    cvmin = min(cvm, na.rm = TRUE)
    idmin = cvm <= cvmin
    lambda.min = max(lambda[idmin], na.rm = TRUE)
    idmin = match(lambda.min, lambda)
    semin = (cvm + cvsd)[idmin]
    id1se = cvm <= semin
    lambda.1se = max(lambda[id1se], na.rm = TRUE)
    id1se = match(lambda.1se, lambda)
    index=matrix(c(idmin,id1se),2,1,dimnames=list(c("min","1se"),"Lambda"))
    list(lambda.min = lambda.min, lambda.1se = lambda.1se, index = index)
}

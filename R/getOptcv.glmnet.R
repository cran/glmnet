getOptcv.glmnet <-
    function (lambda, cvm, cvsd,cvname)
{
    if(match(cvname,c("AUC","C-index"),0))cvm=-cvm
    cvmin = min(cvm, na.rm = TRUE)
    idmin = cvm <= cvmin
    lambda.min = max(lambda[idmin], na.rm = TRUE)
    idmin = match(lambda.min, lambda)
    semin = (cvm + cvsd)[idmin]
    idmin = cvm <= semin
    lambda.1se = max(lambda[idmin], na.rm = TRUE)
    list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

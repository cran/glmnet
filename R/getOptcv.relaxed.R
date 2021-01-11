getOptcv.relaxed <-
    function (statlist,cvname,gamma)
{
    index=matrix(NA,2,2,dimnames=list(c("min","1se"),c("Lambda","Gamma")))
    cvm=lapply(statlist,"[[","cvm")
    nlams=sapply(cvm,length)
    index.lambda=unlist(lapply(nlams,seq))
    cvm=unlist(cvm)
    lambdas=unlist(lapply(statlist,"[[","lambda"))
    cvsd=unlist(lapply(statlist,"[[","cvsd"))
    nzero=unlist(lapply(statlist,"[[","nzero"))
    gammas=rep(gamma,nlams)
    index.gamma = rep(seq(along=gamma),nlams)
    names(lambdas)=NULL
    names(gammas)=NULL
    if(match(cvname,c("AUC","C-index"),0))cvm=-cvm

    cvmin = min(cvm, na.rm = TRUE)
    idmin = cvm <= cvmin
    which=order(lambdas[idmin],gammas[idmin],decreasing=TRUE)[1]
    lambda.min = lambdas[idmin][which]
    gamma.min= gammas[idmin][which]
    index["min","Lambda"]=index.lambda[idmin][which]
    index["min","Gamma"]=index.gamma[idmin][which]
    nzero.min=nzero[idmin][which]
    idmin=seq(along=lambdas)[idmin][which]
    semin = (cvm + cvsd)[idmin]
    idmin = cvm <= semin
    which=order(lambdas[idmin],gammas[idmin],decreasing=TRUE)[1]
    lambda.1se = lambdas[idmin][which]
    gamma.1se= gammas[idmin][which]
    index["1se","Lambda"]=index.lambda[idmin][which]
    index["1se","Gamma"]=index.gamma[idmin][which]
    nzero.1se=nzero[idmin][which]
    list(lambda.min = lambda.min, lambda.1se = lambda.1se, gamma.min=gamma.min, gamma.1se=gamma.1se,nzero.min=nzero.min,nzero.1se=nzero.1se, index=index)
}

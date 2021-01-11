cv.glmnetfit <-function(predmat,y,type.measure,weights,foldid,grouped){
    family=attr(predmat,"family")
    mumat=family$linkinv(predmat)
    nobs=nrow(predmat)# was nrow(mumat), which failed for tweedie instance
    ## initialize from family function. Makes y a vector in case of binomial, and possibly changes weights
    ## Expects nobs to be defined, and creates n and mustart (neither used here)
    ## Some cases expect to see things, so we set it up just to make it work
    etastart=0;mustart=NULL;start=NULL
    eval(family$initialize)
    ##
    ## Just in case this was not done in initialize()
    y <- drop(y)  # we don't like matrix responses
    N = nobs - apply(is.na(predmat), 2, sum)
  cvraw = switch(type.measure,
                 mse = (y - mumat)^2,
                 mae = abs(y - mumat),
                 deviance = family$dev.resids(array(y,dim(mumat)), mumat,1)
                 )
    list(cvraw=cvraw,weights=weights,N=N,type.measure=type.measure,grouped=grouped)
}

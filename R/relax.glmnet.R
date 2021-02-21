#' Fit a relaxed glmnet model
#'
#' @param fit For \code{relax.glmnet} a fitted 'glmnet' object
#' @param maxp a limit on how many relaxed coefficients are allowed. Default is
#' 'n-3', where 'n' is the sample size. This may not be sufficient for
#' non-gaussian familes, in which case users should supply a smaller value.
#' This argument can be supplied directly to 'glmnet'.
#' @param path Since \code{glmnet} does not do stepsize optimization, the Newton
#' algorithm can get stuck and not converge, especially with relaxed fits. With \code{path=TRUE},
#' each relaxed fit on a particular set of variables is computed pathwise using the original sequence
#' of lambda values (with a zero attached to the end). Not needed for Gaussian models, and should not
#' be used unless needed, since will lead to longer compute times. Default is \code{path=FALSE}.
#' appropriate subset of variables
#' @param check.args Should \code{relax.glmnet} make sure that all the data
#' dependent arguments used in creating 'fit' have been resupplied. Default is
#' 'TRUE'.
#' @rdname glmnet
#' @export relax.glmnet
relax.glmnet=function(fit,x,..., maxp = n-3, path=FALSE, check.args=TRUE){

    ## maxp could be a sensitive parameter, especially for glms
    ## This next call fails if you don't supply the crucial arguments
    ## that were used in the creation of "fit". Not needed if called relax.glmnet
    if(check.args)check_dots(fit,...,
                need=c("y","weights","offset","lower.limits","upper.limits"),
                error="used relax.glmnet(),")
    d=dim(x)
    n=d[1];p=d[2]
    ## extract coefs, and find out who is nonzero at each lambda
    a0=fit$a0
    beta=fit$beta
    islistbeta=is.list(beta)
    if(islistbeta){
        ngroups=length(beta)
        nzmat=abs(beta[[1]])
        for(k in 2:ngroups)nzmat=nzmat+abs(beta[[k]])
    }
    else    nzmat=abs(beta)
    nzmat=nzmat>0
    ## Here we go to some trouble to make sure we dont fit the same model more than once
    rgauss=sort(abs(rnorm(p)))
    undex=colSums(rgauss*nzmat)#a unique number for each distinct set of nonzeros
    sundex=sort(unique(undex[undex>0]))
    nzindex=match(sundex,undex)# get a match for each unique
    alldex=match(undex,sundex,0)# for putting them back
    ##
    dev=fit$dev.ratio
    nulldev=fit$nulldev
    excluders=seq(p)
    lam0 = if(path)c(fit$lambda,0)else 0
### Fit a coefficient for each non null model
    internal.parms  <- glmnet.control()
    if (internal.parms$itrace) {
        glmnet.control(itrace = 0) ## disable glmnet reporting
        cat("Relaxed Fits\n")
        pb  <- utils::txtProgressBar(min = 0, max = length(sundex), initial = 0, style = 3)
        on.exit(glmnet.control(itrace = internal.parms$itrace))
        on.exit(close(pb), add = TRUE)
    }
    for(i in seq(along=sundex)){
        ind=nzmat[,nzindex[i]]
        wh=alldex==i
        pi=sum(ind)
        if(pi>maxp)next #can only fit maxp coefficients
        exclude=excluders[!ind]
        ifit=update(fit,x=x,lambda=lam0,exclude=exclude,...,relax=FALSE)
        nlam=ifit$dim[2]
        if (islistbeta) {
            a0[, wh] = ifit$a0[,nlam]
            for (i in 1:ngroups) beta[[i]][, wh] = ifit$beta[[i]][,nlam]
        }
        else {
            a0[wh] = ifit$a0[nlam]
            beta[, wh] = ifit$beta[,nlam]
        }
        dev[wh] = ifit$dev.ratio[nlam]
        if (internal.parms$itrace) utils::setTxtProgressBar(pb, i)
    }
    if (internal.parms$itrace) {
        utils::setTxtProgressBar(pb = pb, value = length(sundex))
    }

    fitg=fit
    omit=colSums(nzmat)>maxp
    if(any(omit)){
        if(islistbeta){
            a0=a0[,!omit]
            for(i in 1:ngroups)beta[[i]]=beta[[i]][,!omit]
        }
        else{
            beta=beta[,!omit]
            a0=a0[!omit]
            }
        dev=dev[!omit]
        fit$lambda=fit$lambda[!omit]
        fit$df=fit$df[!omit]
        }
    fit$beta=beta
    fit$a0=a0
    fit$dev.ratio=dev
    class(fitg)=c("relaxed",class(fitg))
    fitg$relaxed=fit
    fitg
    }





blend.relaxed=function(fit,gamma,extend=TRUE,epsilon=0.00001){
###
###gamma must be a single value
    if(gamma==1){
        class(fit)=class(fit)[-1] #drop the relaxed status
        fit$relaxed=NULL
        return(fit)
    }
    gamma=max(gamma,epsilon)
### We do the above so that even of gamma=0, and extend=TRUE, we fill in the coefficients
### First see which if any of glmnet coefficients are missing on relaxed
    which=match(fit$lambda,fit$relaxed$lambda,0)
    if(extend)which[which==0]=max(which)
    beta=fit$beta
    a0=fit$a0
    fitr=fit$relaxed
    betar=fitr$beta
    a0r=fitr$a0
    islistbeta=is.list(beta)
    if(!islistbeta){
        fitr$beta=gamma*beta[,which>0] +(1-gamma)*betar[,which]
        fitr$a0=gamma*a0[which>0]+(1-gamma)*a0r[which]
    }
    else{
        ngroups=length(beta)
        for(i in 1:ngroups)
            beta[[i]]=gamma*beta[[i]][,which>0]+(1-gamma)*betar[[i]][,which]
        fitr$beta=beta
        fitr$a0=gamma*a0[,which>0]+(1-gamma)*a0r[,which]
    }
    fitr$dev.ratio=gamma*fit$dev.ratio[which>0]+(1-gamma)*fitr$dev.ratio[which]
    fitr$df=fit$df[which>0]
    fitr$lambda=fit$lambda[which>0]
    fitr
    }

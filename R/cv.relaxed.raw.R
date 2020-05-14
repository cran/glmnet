cv.relaxed.raw <-
    function (x, y, weights, offset, lambda, type.measure, nfolds, foldid, alignment,grouped, keep,
              parallel, trace.it, glmnet.call, cv.call, gamma, ...)
{

###This next function may be redundant
  relaxglmnet=function(x,...){
    fit=glmnet(x,...)
    relax.glmnet(fit,x,...,check.args=FALSE)
     }
    gamma=sort(checkgamma.relax(gamma))
### The next bit is to make sure we always include lasso in CV, as well as the gamma series
    lengamma=length(gamma)
    if(max(gamma)<1.0){
        gamma=c(gamma,1.0)
        }

    glmnet.call$relax=TRUE
  if (trace.it) cat("Training\n")
  glmnet.object = relaxglmnet(x=x, y=y, weights = weights, offset = offset,
    lambda = lambda, trace.it=trace.it,...)
  glmnet.object$call = glmnet.call
  subclass=class(glmnet.object)[[2]]# it is of class c("relaxed","subtype","glmnet")
  type.measure=cvtype(type.measure,subclass)
  is.offset = glmnet.object$offset
###Next line is commented out so each call generates its own lambda sequence
# lambda=glmnet.object$lambda
 if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
    nz = predict(glmnet.object, type = "nonzero")
    nz = sapply(nz, function(x) sapply(x, length))
    nz = ceiling(apply(nz, 1, median))
  }
  else nz = sapply(predict(glmnet.object, type = "nonzero"),
         length)

  outlist = as.list(seq(nfolds))
  N=nrow(x)
    if (parallel) {
#  if (parallel && require(foreach)) {
    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar%
    {
      which = foldid == i
#      if (is.matrix(y))
      if (length(dim(y))>1)
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset)
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      relaxglmnet(x=x[!which, , drop = FALSE], y=y_sub, lambda = lambda,
             offset = offset_sub, weights = weights[!which],
             ...)
    }
  }
  else {
      for (i in seq(nfolds)) {
      if (trace.it) cat(sprintf("Fold: %d/%d\n", i, nfolds))
      which = foldid == i
      if (is.matrix(y))
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset)
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      outlist[[i]] = relaxglmnet(x=x[!which, , drop = FALSE],
               y=y_sub, lambda = lambda, offset = offset_sub,
               weights = weights[!which],trace.it=trace.it, ...)
    }
  }
  lambda = glmnet.object$lambda
  class(outlist)=paste0(subclass,"list")
  predmatlist=as.list(gamma)
  names(predmatlist)=paste("g",gamma,sep=":")
  outstuff=cvstufflist=predmatlist
### Even though the following is innefficient, It makes the code more modular
###  fun = paste("cv", subclass, sep = ".")
  for(i in seq(along=gamma)){
      predob=buildPredmat(outlist,
                         lambda=lambda,x=x,offset=offset,foldid=foldid,alignment=alignment,
                         y=y,weights=weights,grouped=grouped,
                         gamma=gamma[i],family=family(glmnet.object)
                         )
#      if(subclass=="glmnetfit")attr{predob,"family")=glmnet.object$family #special case
      predmatlist[[i]]=predob
      fun = paste("cv", subclass, sep = ".")
      cvstufflist[[i]] = do.call(fun, list(predmatlist[[i]],y,type.measure,weights,foldid,grouped))
  }
  grouped=all(sapply(cvstufflist,"[[","grouped"))
    if ((N/nfolds < 3) && grouped) {
    warning("Option grouped=FALSE enforced in cv.relaxed, since < 3 observations per fold",
            call. = FALSE)
    grouped = FALSE
    }
  for(i in seq(along=gamma)) outstuff[[i]]=cvstats(cvstufflist[[i]],foldid,nfolds,lambda,nz,grouped)
  cvn = cvstufflist[[1]]$type.measure
    cvname=names(cvn);names(cvname)=cvn# to be compatible with earlier version; silly, I know
    out=outstuff[[length(gamma)]] # lasso stats
    out = c(out,list(call=cv.call,name = cvname, glmnet.fit = glmnet.object))
    if(lengamma<length(gamma)){
        outstuff=outstuff[1:lengamma]
        gamma=gamma[1:lengamma]
        }
    relaxed=list(statlist=outstuff,gamma=gamma)
    lamin=getOptcv.relaxed(outstuff,cvname,gamma)
    relaxed=  c(relaxed, as.list(lamin))# for relaxed part
    lamin=with(out,getOptcv.glmnet(lambda, cvm, cvsd, cvname))
    out = c(out, as.list(lamin)) # for lasso part
  if (keep)
      out = c(out, list(fit.preval = predmatlist, foldid = foldid))
    out$relaxed=relaxed
  class(out) = c("cv.relaxed","cv.glmnet")
  out
}


cv.glmnet.raw <-
  function (x, y, weights, offset, lambda, type.measure, nfolds, foldid, alignment,grouped, keep,
              parallel, trace.it, glmnet.call, cv.call, ...)
{


  if (trace.it) cat("Training\n")
  glmnet.object = glmnet(x, y, weights = weights, offset = offset,
    lambda = lambda, trace.it=trace.it,...)
  glmnet.object$call = glmnet.call
  subclass=class(glmnet.object)[[1]]
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
      glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda,
             offset = offset_sub, weights = weights[!which],
             ...)
    }
  }
  else {
      for (i in seq(nfolds)) {
          if (trace.it) cat(sprintf("Fold: %d/%d\n", i, nfolds))
        which = foldid == i
        if (length(dim(y))>1)
            y_sub = y[!which, ]
        else y_sub = y[!which]
      if (is.offset)
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      outlist[[i]] = glmnet(x[!which, , drop = FALSE],
               y_sub, lambda = lambda, offset = offset_sub,
               weights = weights[!which],trace.it=trace.it, ...)
    }
  }
  lambda = glmnet.object$lambda
  class(outlist)=paste0(subclass,"list")
  predmat=buildPredmat(outlist,lambda,x,offset,foldid,alignment,y=y,weights=weights,
                       grouped=grouped,type.measure=type.measure,family=family(glmnet.object))
### we include type.measure for the special case of coxnet with the deviance vs C-index discrepancy
### family is included for the new GLM crowd
### Next we compute the measures
#    if(subclass=="glmnetfit") attr(predmat,"family")=glmnet.object$family
    fun = paste("cv", subclass, sep = ".")
    cvstuff = do.call(fun, list(predmat,y,type.measure,weights,foldid,grouped))

  grouped=cvstuff$grouped
    if ((N/nfolds < 3) && grouped) {
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",
            call. = FALSE)
    grouped = FALSE
    }

  out=cvstats(cvstuff,foldid,nfolds,lambda,nz,grouped)
  cvname = names(cvstuff$type.measure)
  names(cvname)=cvstuff$type.measure# to be compatible with earlier version; silly, I know
  out = c(out,list(call=cv.call,name = cvname, glmnet.fit = glmnet.object))
  if (keep)
    out = c(out, list(fit.preval = predmat, foldid = foldid))
  lamin=with(out,getOptcv.glmnet(lambda, cvm, cvsd, cvname))
  obj = c(out, as.list(lamin))
  class(obj) = "cv.glmnet"
  obj
}


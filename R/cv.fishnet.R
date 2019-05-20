cv.fishnet <-
  function (outlist, lambda, x, y, weights, offset, foldid, type.measure,
            grouped, keep = FALSE,alignment="lambda")
{
  if (!is.null(offset)) {
    is.offset = TRUE
    offset = drop(offset)
  }
  else is.offset = FALSE

  devi = function(y, eta) {
    deveta = y * eta - exp(eta)
    devy = y * log(y) - y
    devy[y == 0] = 0
    2 * (devy - deveta)
  }
  predmat = matrix(NA, length(y), length(lambda))
  nfolds = max(foldid)
  nlams = double(nfolds)
  nlambda=length(lambda)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    if (is.offset)
      off_sub = offset[which]
    preds = switch(alignment,
                   fraction=predict(fitobj, x[which, , drop = FALSE], newoffset = off_sub),
                   lambda=predict(fitobj, x[which, , drop = FALSE], s=lambda, newoffset = off_sub)
                   )
    nlami = min(ncol(preds),nlambda)
    predmat[which, seq(nlami)] = preds[,seq(nlami)]
    if(nlami<nlambda)predmat[which,seq(from=nlami,to=nlambda)]=preds[,nlami]
    nlams[i] = nlambda
  }
  N = length(y) - apply(is.na(predmat), 2, sum)
  cvraw = switch(type.measure,
                 mse = (y - exp(predmat))^2,
                 mae = abs(y - exp(predmat)),
                 deviance = devi(y, predmat)
                 )
  if ((length(y)/nfolds < 3) && grouped) {
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",
            call. = FALSE)
    grouped = FALSE
  }
  if (grouped) {
    cvob = cvcompute(cvraw, weights, foldid, nlams)
    cvraw = cvob$cvraw
    weights = cvob$weights
    N = cvob$N
  }
  cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
  cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
    w = weights, na.rm = TRUE)/(N - 1))
  out = list(cvm = cvm, cvsd = cvsd, type.measure=type.measure)
  if (keep)
    out$fit.preval = predmat
  out
}

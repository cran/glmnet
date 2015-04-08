cv.multnet <-
  function (outlist, lambda, x, y, weights, offset, foldid, type.measure, 
            grouped, keep = FALSE) 
{
  typenames = c(mse = "Mean-Squared Error", mae = "Mean Absolute Error", 
    deviance = "Multinomial Deviance", class = "Misclassification Error")
  if (type.measure == "default") 
    type.measure = "deviance"
  if (!match(type.measure, c("mse", "mae", "deviance", "class"), 
             FALSE)) {
    warning("Only 'deviance', 'class',  'mse' or 'mae'  available for multinomial models; 'deviance' used")
    type.measure = "deviance"
  }
  prob_min = 1e-05
  prob_max = 1 - prob_min
  nc = dim(y)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
  }
  else nc = nc[2]
  if (!is.null(offset)) 
    is.offset = TRUE
  else is.offset = FALSE

  mlami=max(sapply(outlist,function(obj)min(obj$lambda)))
  which_lam=lambda >= mlami

  predmat = array(NA, c(nrow(y), nc, length(lambda)))
  nfolds = max(foldid)
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    if (is.offset) 
      off_sub = offset[which, , drop = FALSE]
    preds = predict(fitobj, x[which, , drop = FALSE], s=lambda[which_lam],offset = off_sub, 
      type = "response")
    nlami = sum(which_lam)
    predmat[which, , seq(nlami)] = preds
    nlams[i] = nlami
  }
  ywt = apply(y, 1, sum)
  y = y/ywt
  weights = weights * ywt
  N = nrow(y) - apply(is.na(predmat[, 1, ]), 2, sum)
  bigY = array(y, dim(predmat))
  cvraw = switch(type.measure, mse = apply((bigY - predmat)^2, 
                                 c(1, 3), sum), mae = apply(abs(bigY - predmat), c(1, 
                                                  3), sum), deviance = {
                                                    predmat = pmin(pmax(predmat, prob_min), prob_max)
                                                    lp = bigY * log(predmat)
                                                    ly = bigY * log(bigY)
                                                    ly[bigY == 0] = 0
                                                    apply(2 * (ly - lp), c(1, 3), sum)
                                                  }, class = {
                                                    classid = as.numeric(apply(predmat, 3, glmnet_softmax))
                                                    yperm = matrix(aperm(bigY, c(1, 3, 2)), ncol = nc)
                                                    matrix(1 - yperm[cbind(seq(classid), classid)], ncol = length(lambda))
                                                  })
  if ((nrow(y)/nfolds < 3) && grouped) {
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
  out = list(cvm = cvm, cvsd = cvsd, name = typenames[type.measure])
  if (keep) 
    out$fit.preval = predmat
  out
}

predict.multnet <-
  function (object, newx, s = NULL, type = c("link", "response",
                                      "coefficients", "class", "nonzero"), exact = FALSE, newoffset,
            ...)
{
  type = match.arg(type)
    ###multnet  is very different, so we treat it separately
  if(missing(newx)){
    if(!match(type,c("coefficients","nonzero"),FALSE))stop("You need to supply a value for 'newx'")
  }
  if(exact&&(!is.null(s))){
    lambda=object$lambda
    which=match(s,lambda,FALSE)
    if(!all(which>0)){
      lambda=unique(rev(sort(c(s,lambda))))
      check_dots(object,...)# This fails if you don't supply the crucial arguments
      object=update(object,lambda=lambda,...)
    }
  }

  a0 = object$a0
  rownames(a0) = rep("(Intercept)", nrow(a0))
  nbeta = object$beta
  klam = dim(a0)
  nclass = klam[[1]]
  nlambda = length(s)
  if (!is.null(s)) {
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    for (i in seq(nclass)) {
      kbeta = methods::rbind2(a0[i, , drop = FALSE], nbeta[[i]])#was rbind2
      vnames = dimnames(kbeta)[[1]]
      dimnames(kbeta) = list(NULL, NULL)
      kbeta = kbeta[, lamlist$left, drop = FALSE] %*% Diagonal(x=lamlist$frac) +
        kbeta[, lamlist$right, drop = FALSE] %*% Diagonal(x=1 - lamlist$frac)
      dimnames(kbeta) = list(vnames, paste(seq(along = s)))
      nbeta[[i]] = kbeta
    }
  }
  else {
    for (i in seq(nclass)) nbeta[[i]] = methods::rbind2(a0[i, ], nbeta[[i]])#was rbind2
    nlambda=length(object$lambda)
  }
  if (type == "coefficients")
    return(nbeta)
  if (type == "nonzero")
    if(object$grouped)return(nonzeroCoef(object$beta[[1]],bystep=TRUE))
    else
      return(lapply(nbeta, function(x) nonzeroCoef(x[-1, ,
                                                   drop = FALSE], bystep = TRUE)))
  dd = dim(newx)
   if (inherits(newx, "sparseMatrix"))
    newx = as(newx, "dgCMatrix")
  npred = dd[[1]]
  dn = list(names(nbeta), dimnames(nbeta[[1]])[[2]], dimnames(newx)[[1]])
   dp = array(0, c(nclass, nlambda, npred), dimnames = dn)
  for (i in seq(nclass)) {
    fitk = cbind2(1, newx) %*% (nbeta[[i]])
    dp[i, , ] = dp[i, , ] + t(as.matrix(fitk))
  }
  if (object$offset) {
    if (missing(newoffset))
      stop("No newoffset provided for prediction, yet offset used in fit of glmnet",
           call. = FALSE)
    if (!is.matrix(newoffset) || dim(newoffset)[[2]] != nclass)
      stop(paste("Dimension of newoffset should be", npred, "x",
                 nclass), call. = FALSE)
    toff = t(newoffset)
    for (i in seq(nlambda)) {
      dp[, i, ] = dp[, i, ] + toff
    }
  }
  switch(type, response = {
    pp = exp(dp)
    psum = apply(pp, c(2, 3), sum)
    aperm(pp/rep(psum, rep(nclass, nlambda * npred)), c(3,
                                                        1, 2))
  }, link = aperm(dp, c(3, 1, 2)), class = {
    dp = aperm(dp, c(3, 1, 2))
    apply(dp, 3, glmnet_softmax)
  })
}

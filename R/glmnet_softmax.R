glmnet_softmax <-
  function (x) 
{
  d <- dim(x)
  nas=apply(is.na(x),1,any)
  if(any(nas)){
    pclass=rep(NA,d[1])
    if(sum(nas)<d[1]){
      pclass2=glmnet_softmax(x[!nas,])
      pclass[!nas]=pclass2
      if(is.factor(pclass2))pclass=factor(pclass,levels = seq(d[2]),labels=levels(pclass2))
    }
  }
  else{
    maxdist <- x[, 1]
    pclass <- rep(1, d[1])
    for (i in seq(2, d[2])) {
      l <- x[, i] > maxdist
      pclass[l] <- i
      maxdist[l] <- x[l, i]
    }
    dd <- dimnames(x)[[2]]
    pclass <- if (is.null(dd) || !length(dd)) 
      pclass
    else factor(pclass, levels = seq(d[2]), labels = dd)
  }
  pclass
}


glmnet_softmax <-
  function (x,ignore_labels=FALSE)
{
  d <- dim(x)
  dd <- dimnames(x)[[2]]
  if(is.null(dd) || !length(dd)) ignore_labels=TRUE

  nas=apply(is.na(x),1,any)
  if(any(nas)){
    pclass=rep(NA,d[1])
    if(sum(nas)<d[1]){
      pclass2=glmnet_softmax(x[!nas,],ignore_labels)
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
    if(!ignore_labels) pclass=factor(pclass, levels = seq(d[2]), labels = dd)
  }
  pclass
}


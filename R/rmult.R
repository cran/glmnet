rmult <-
  function(p){
    x=t(apply(p,1,function(x)rmultinom(1,1,x)))
    x=x%*%seq(ncol(p))
    drop(x)
  }


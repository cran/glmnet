nonzeroCoef = function (beta, bystep = FALSE) 
{
### bystep = FALSE means which variables were ever nonzero
### bystep = TRUE means which variables are nonzero for each step
  nr=nrow(beta)
  if (nr == 1) {#degenerate case
    if (bystep) 
      apply(beta, 2, function(x) if (abs(x) > 0) 
            1
      else NULL)
    else {
      if (any(abs(beta) > 0)) 
        1
      else NULL
    }
  }
  else {
    beta=abs(beta)>0 # this is sparse
    which=seq(nr)
    ones=rep(1,ncol(beta))
    nz=as.vector((beta%*%ones)>0)
    which=which[nz]
    if (bystep) {
      beta=as.matrix(beta[which,])
      nzel = function(x, which) if (any(x)) 
        which[x]
      else NULL
      apply(beta, 2, nzel, which)
    }
    else which
  }
}

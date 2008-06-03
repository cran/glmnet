lambda.interp=function(lambda,s){
### lambda is the index sequence that is produced by the model
### s is the new vector at which evaluations are required.
### the value is a vector of left and right indices, and a vector of fractions.
### the new values are interpolated bewteen the two using the fraction
### Note: lambda decreases. you take:
### sfrac*left+(1-sfrac*right)
  
  if(length(lambda)==1){# degenerate case of only one lambda
    nums=length(s)
    left=rep(1,nums)
    right=left
    sfrac=rep(1,nums)
  }
  else{
      s[s > max(lambda)] = max(lambda)
      s[s < min(lambda)] = min(lambda)
      k=length(lambda)
      sfrac <- (lambda[1]-s)/(lambda[1] - lambda[k])
      lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
      coord <- approx(lambda, seq(lambda), sfrac)$y
      left <- floor(coord)
      right <- ceiling(coord)
      sfrac=(sfrac-lambda[right])/(lambda[left] - lambda[right])
      sfrac[left==right]=1
    }
list(left=left,right=right,frac=sfrac)
}

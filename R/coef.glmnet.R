coef.glmnet=function(object,s=object$lambda,exact=FALSE,...)
  predict(object,s=s,type="coefficients",exact=exact)

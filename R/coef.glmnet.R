coef.glmnet=function(object,s=NULL,exact=FALSE,...)
  predict(object,s=s,type="coefficients",exact=exact)

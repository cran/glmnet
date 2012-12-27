coef.glmnet=function(object,s=NULL,exact=TRUE,...)
  predict(object,s=s,type="coefficients",exact=exact)

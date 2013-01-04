predict.lognet=function(object,newx,s=NULL,type=c("link","response","coefficients","class","nonzero"),exact=FALSE,offset,...){
  type=match.arg(type)
    ### remember that although the fortran lognet makes predictions
  ### for the first class, we make predictions for the second class
  ### to avoid confusion with 0/1 responses.
  ### glmnet flipped the signs of the coefficients 
nfit=NextMethod("predict")
  switch(type,
         response={
           pp=exp(-nfit)
           1/(1+pp)
         },
         class={
           cnum=ifelse(nfit>0,2,1)
           clet=object$classnames[cnum]
           if(is.matrix(cnum))clet=array(clet,dim(cnum),dimnames(cnum))
           clet
         },
         nfit
         )
}  
    

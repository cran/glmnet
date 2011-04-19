predict.fishnet=function(object,newx,s=NULL,type=c("link","response","coefficients","nonzero"),exact=FALSE,offset,...){
  type=match.arg(type)
nfit=NextMethod("predict")
  switch(type,
         response=exp(nfit),
         nfit
         )
}  

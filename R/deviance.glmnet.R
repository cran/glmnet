deviance.glmnet=function(object,...){
dev=object$dev
nulldev=object$nulldev
(1-dev)*nulldev
}

#' compute C index for a Cox model
#'
#' Computes Harrel's C index for predictions from a \code{"coxnet"} object.
#'
#' Computes the concordance index, taking into account censoring.
#'
#' @param pred Predictions from a \code{"coxnet"} object
#' @param y a survival response object - a matrix with two columns "time" and
#' "status"; see documentation for "glmnet"
#' @param weights optional observation weights
#' @author Rob Tibshirani\cr Maintainer: Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{cv.glmnet}
#' @references Harrel Jr, F. E. and Lee, K. L. and Mark, D. B. (1996)
#' \emph{Tutorial in biostatistics: multivariable prognostic models: issues in
#' developing models, evaluating assumptions and adequacy, and measuring and
#' reducing error}, Statistics in Medicine, 15, pages 361--387.
#' @keywords Cox models cross-validation
#' @examples
#'
#' set.seed(10101)
#' N = 1000
#' p = 30
#' nzc = p/3
#' x = matrix(rnorm(N * p), N, p)
#' beta = rnorm(nzc)
#' fx = x[, seq(nzc)] %*% beta/3
#' hx = exp(fx)
#' ty = rexp(N, hx)
#' tcens = rbinom(n = N, prob = 0.3, size = 1)  # censoring indicator
#' y = cbind(time = ty, status = 1 - tcens)  # y=Surv(ty,1-tcens) with library(survival)
#' fit = glmnet(x, y, family = "cox")
#' pred = predict(fit, newx = x)
#' Cindex(pred, y)
#' cv.glmnet(x, y, family = "cox", type.measure = "C")
#'
#' @export Cindex
Cindex=function(pred,y,weights=rep(1,nrow(y))){
###  Written by Rob Tibshirani; edited by Trevor Hastie
###  This function gives identical results to the summary function from coxph,
###  and the concordance.index function in survcomp (with their default settings),
###  and no ties in yhat.
###  Works with ties in y. But does not agree with latter when yhat has ties.
###  There are conflicting definitions for c-index in this case
###
###  Formula used  is
###  Concordance = (num. all concordant pairs + num. tied pairs/2)/(num. total pairs including ties).
###  with weights,  weights used are average wts for each pair (??)

    status=y[,"status"]
    ty=y[,"time"]
    risksets=which(status==1)

    weights=length(weights)*weights/sum(weights)
    fun=function(riskset,ty,pred,weights){
           total=concordant=0
        i=riskset
        rest=which(ty>ty[i])
        if(length(rest)>0){

            ww=(weights[rest]+weights[i])/2
             total=sum(ww)
          concordant = 1.0*sum(ww*(pred[rest]<pred[i]))+0.5*sum(ww*(pred[rest]==pred[i]))
        }

    return(c(concordant,total))
}

  out=  sapply(risksets,fun,ty,pred,weights)
    cindex=sum(out[1,])/sum(out[2,])

    return(cindex)
}

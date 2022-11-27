#' fit a glm with all the options in \code{glmnet}
#'
#' Fit a generalized linear model as in \code{glmnet} but unpenalized. This
#' allows all the features of \code{glmnet} such as sparse x, bounds on
#' coefficients, offsets, and so on.
#'
#' This is essentially the same as fitting a "glmnet" model with a single value
#' \code{lambda=0}, but it avoids some edge cases. CAVEAT: If the user tries a
#' problem with N smaller than or close to p for some models, it is likely to
#' fail (and maybe not gracefully!) If so, use the \code{path=TRUE} argument.
#'
#' @param x input matrix
#' @param ... Most other arguments to glmnet that make sense
#' @param path  Since \code{glmnet} does not do stepsize optimization, the Newton
#' algorithm can get stuck and not converge, especially with unpenalized fits. With \code{path=TRUE},
#' the  fit computed with pathwise lasso regularization. The current implementation does this twice:
#' the first time to get the lambda sequence, and the second time with a zero attached to the end).
#' Default is \code{path=FALSE}.
#' @return It returns an object of class "bigGlm" that inherits from class
#' "glmnet". That means it can be predicted from, coefficients extracted via
#' \code{coef}. It has its own print method.
#' @author Trevor Hastie\cr Maintainer: Trevor Hastie
#' \email{hastie@@stanford.edu}
#' @seealso \code{print}, \code{predict}, and \code{coef} methods.
#' @keywords models regression
#' @examples
#'
#' # Gaussian
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit1 = bigGlm(x, y)
#' print(fit1)
#'
#' fit2=bigGlm(x,y>0,family="binomial")
#' print(fit2)
#' fit2p=bigGlm(x,y>0,family="binomial",path=TRUE)
#' print(fit2p)
#'
#' @export bigGlm
bigGlm=function (x,..., path=FALSE){
### fit a model using x, y  and the other glmnet specs in ... unpenalized
    thiscall = match.call(glmnet)
### Next few lines simply in case x is a vector, in which case glmnet complains
    bx=cbind(x,1)
    d=dim(bx)
    n=d[1]
    p=d[2]-1
    bx[1,p+1]=0 # so its not a constant
    arglist = list(...)
    arglist$x = bx
    arglist$exclude=p+1
    fixbeta <- function(beta,p,nlam=1){
        beta <- beta[1:p, nlam , drop = FALSE]
        dn <- dimnames(beta)
        dn <- list(dn[[1]],NULL)
        dimnames(beta) <- dn
        return(beta)
        }
    if(!path){
        arglist$lambda=0
        fit = do.call("glmnet", arglist)
        if (is.list(fit$beta)){
            fit$beta = lapply(fit$beta, fixbeta, p=p)
            dimnames(fit$a0) <- NULL
            }
        else {
                fit$beta = fixbeta(fit$beta, p=p)
                names(fit$a0) <- NULL
                }
    }
    else
    {
        fit = do.call("glmnet", arglist)# get the lambda sequence
        lam0=c(fit$lambda,0)
        arglist$lambda=lam0
        fit = do.call("glmnet", arglist)
        nlam=length(fit$lambda)#may have stopped early
        if (is.list(fit$beta)){
            fit$a0=fit$a0[,nlam,drop=FALSE]
            dimnames(fit$a0) <- NULL
            fit$beta = lapply(fit$beta, fixbeta, p=p, nlam=nlam)
            dfmat=fit$dfmat
            if(!is.null(dfmat))fit$dfmat=dfmat[,nlam,drop=FALSE]
        }
        else {
            fit$a0=fit$a0[nlam]
            names(fit$a0) <- NULL
            fit$beta = fixbeta(fit$beta, p=p, nlam=nlam)
            }
        fit$df=p
        fit$lambda=0
        fit$dev.ratio=fit$dev.ratio[nlam]
    }
    fit$dim=c(p,1)
    fit$call = thiscall
    class(fit) = c("bigGlm", class(fit))
    fit
}



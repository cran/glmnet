#' @method plot relaxed
#' @param gamma Value of the mixing parameter for a "relaxed" fit
#' @rdname plot.glmnet
#' @export
plot.relaxed <-
    function(x,xvar = c("lambda", "dev"), label = FALSE, gamma=1,
             ...) {
        xvar=match.arg(xvar)
        if(any(wh<-gamma<0)){
            warning("negative gamma values ignored")
            gamma=gamma[!wh]
        }
        if(any(wh<-gamma>1)){
            warning("gamma values larger than 1 ignored")
            gamma=gamma[!wh]
        }
        if(!length(gamma))stop("no valid values of gamma")
        for(i in seq(along=length(gamma))){
            bfit=blend.relaxed(x,gamma=gamma[i],extend=FALSE)
            plot(bfit,xvar=xvar,label=label,...)
            }
    }


#' @method print relaxed
#' @export
print.relaxed <-
function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n")
    cat("Relaxed\n\n")
    which=match(x$relaxed$lambda,x$lambda,0)
    rdev=rep(NA,length(x$lambda))
    rdev[which]=round(x$relaxed$dev.ratio*100, 2)
    out=data.frame(Df = x$df, `%Dev` = round(x$dev.ratio*100, 2), `%Dev R`=rdev,
                    Lambda = signif(x$lambda, digits),check.names=FALSE,row.names=seq(along=rdev))
    class(out)=c("anova",class(out))
    print(out)
}


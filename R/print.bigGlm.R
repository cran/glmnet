#' @method print bigGlm
#' @export
print.bigGlm <-
function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
     out=data.frame(Df = x$df, `%Dev` = round(x$dev.ratio*100, 2),
                    Lambda = 0,check.names=FALSE,row.names=seq(along=x$df))
class(out)=c("anova",class(out))
     print(out)
}


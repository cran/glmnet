#' @method print cv.relaxed
#' @export
print.cv.relaxed <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    cat("Measure:", x$name, "\n\n")
    x = x$relaxed
    optlams = c(x$lambda.min, x$lambda.1se)
    wg1 = match(x$gamma.min, x$gamma)
    wl1 = match(x$lambda.min, x$statlist[[wg1]]$lambda)
    s1 = with(x$statlist[[wg1]], c(x$gamma.min,wg1, x$lambda.min,wl1,
        cvm[wl1], cvsd[wl1], x$nzero.min))
    wg2 = match(x$gamma.1se, x$gamma)
    wl2 = match(x$lambda.1se, x$statlist[[wg2]]$lambda)
    s2 = with(x$statlist[[wg2]], c(x$gamma.1se,wg2, x$lambda.1se,wl2,
        cvm[wl2], cvsd[wl2], x$nzero.1se))
    mat = rbind(s1, s2)
    dimnames(mat) = list(c("min", "1se"), c("Gamma","Index", "Lambda","Index",
        "Measure", "SE", "Nonzero"))
    mat = data.frame(mat, check.names = FALSE)
    class(mat) = c("anova", class(mat))
    print(mat, digits = digits)
}

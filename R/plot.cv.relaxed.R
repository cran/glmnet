#' method plot cv.relaxed
#' @importFrom shape colorlegend
#' @param se.bands Should shading be produced to show standard-error bands;
#' default is \code{TRUE}
#' @rdname plot.cv.glmnet
#' @export
plot.cv.relaxed <- function (x, se.bands=TRUE,...)
{
    xr=x$relaxed
    oldpar = par(mar = c(4, 4, 3, 4))
    on.exit(par(oldpar))
    statlist=xr$statlist
    gamma=xr$gamma
    ngamma=length(gamma)
    ylim=range(unlist(lapply(statlist,"[[","cvm")))
    if(se.bands){
        cvup=lapply(statlist,"[[","cvup")
        cvlo=lapply(statlist,"[[","cvlo")
        ylim=range(ylim,unlist(cvup),unlist(cvlo))
        }
    xlim=log(range(unlist(lapply(statlist,"[[","lambda")))+0.00001)
    cvcolors = rainbow(ngamma, start = .1, end = 1)
     with(statlist[[ngamma]],plot(log(lambda), cvm, type = "n", xlab = expression(Log(lambda)),
                                 ylab = x$name,ylim=ylim,xlim=xlim))
    if(se.bands){
        for (i in seq(ngamma))
            with(statlist[[i]],polygon(c(log(lambda),rev(log(lambda))),c(cvup,rev(cvlo)),
                                       col="floralwhite",border="antiquewhite"))
    }
    for (i in seq(ngamma))
        with(statlist[[i]],lines(log(lambda), cvm, lwd = 2, col = cvcolors[i]))
    mins=log(c(xr$lambda.min,xr$lambda.1se))
    abline(v=mins,lty=3)
    dof= statlist[[1]]$nzero
    lambda=statlist[[1]]$lambda
    axis(side = 3, at = log(lambda), labels = paste(dof),
         tick = FALSE, line = 0)

#    dof=c(x$nzero.min,x$nzero.1se)
#    axis(3, at = mins, labels = format(dof), tick=FALSE, line = 0,cex.axis=.9)
    colorlegend(posy = c(0.2, 0.8), posx = c(0.93, 0.945)-.03, col = rainbow(ngamma,
        start = 0.1, end = 1), zlim = c(0, 1), zval = gamma , main = expression(gamma),digit=2)
    invisible()
}

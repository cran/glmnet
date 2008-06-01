`softmax` <-
function (x, gap = FALSE) 
{
    d <- dim(x)
    maxdist <- x[, 1]
    pclass <- rep(1, d[1])
    for (i in seq(2, d[2])) {
        l <- x[, i] > maxdist
        pclass[l] <- i
        maxdist[l] <- x[l, i]
    }
    dd <- dimnames(x)[[2]]
    if (gap) {
        x <- abs(maxdist - x)
        x[cbind(seq(d[1]), pclass)] <- drop(x %*% rep(1, d[2]))
        gaps <- do.call("pmin", data.frame(x))
    }
    pclass <- if (is.null(dd) || !length(dd)) 
        pclass
    else factor(pclass, levels = seq(d[2]), labels = dd)
    if (gap) 
        list(class = pclass, gaps = gaps)
    else pclass
}

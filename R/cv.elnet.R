cv.elnet <-function(predmat,y,type.measure,weights,foldid,grouped){

    N = length(y) - apply(is.na(predmat), 2, sum)
    cvraw = switch(type.measure,
                   mse = (y - predmat)^2,
                   deviance = (y - predmat)^2,
                   mae = abs(y - predmat))
    list(cvraw=cvraw,weights=weights,N=N,type.measure=type.measure,grouped=grouped)

}

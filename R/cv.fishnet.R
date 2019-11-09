cv.fishnet <-function(predmat,y,type.measure,weights,foldid,grouped){

    devi = function(y, eta) {
    deveta = y * eta - exp(eta)
    devy = y * log(y) - y
    devy[y == 0] = 0
    2 * (devy - deveta)
      }

  N = length(y) - apply(is.na(predmat), 2, sum)
  cvraw = switch(type.measure,
                 mse = (y - exp(predmat))^2,
                 mae = abs(y - exp(predmat)),
                 deviance = devi(y, predmat)
                 )
    list(cvraw=cvraw,weights=weights,N=N,type.measure=type.measure,grouped=grouped)
}

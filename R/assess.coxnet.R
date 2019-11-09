assess.coxnet <- function(predmat,y,type.measure,weights,foldid,grouped){
### This is a hack, because CV for coxnet is complex
if(type.measure=="C")return(cv.coxnet(predmat,y,type.measure,weights,foldid,grouped))
    nlambda=ncol(predmat)
    N=nrow(predmat)
    cvraw = matrix(NA, 1, nlambda)
    for (j in seq(nlambda)) {
        cvraw[1, j] = coxnet.deviance(predmat[,j],y, weights=weights)
            }
    list(cvraw=cvraw,weights=N,N=N,type.measure=type.measure,grouped=FALSE)
}

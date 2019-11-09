#' @method print confusion.table
#' @export
print.confusion.table=function(x,digits = max(3, getOption("digits") - 3), ...){

    ndn=names(dimnames(x))
    rowtot=rowSums(x)
    coltot=colSums(x)
    tot=sum(coltot)
    ncorrect=sum(diag(x))
    correct=(ncorrect)/tot

    x=cbind(x,Total=rowtot)
    x=rbind(x,Total=c(coltot,tot))
    dn=dimnames(x)
    names(dn)=ndn
    dimnames(x)=dn
    print(x,digits=digits,...)
    cat("\n Percent Correct: ",format(round(correct,digits)),"\n")
    invisible()
}

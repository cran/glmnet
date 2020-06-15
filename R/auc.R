auc=function(y,prob,w){
    if(missing(w))
        concordance(y~prob)$concordance
    else concordance(y~prob,weights=w)$concordance
    }

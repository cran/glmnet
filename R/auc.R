auc=function(y,prob,w){
    if(missing(w)){
        rprob=rank(prob)
        n1=sum(y);n0=length(y)-n1
        u=sum(rprob[y==1])-n1*(n1+1)/2
        exp(log(u) - log(n1) - log(n0))
    }
    else{
        rprob=runif(length(prob))
        op=order(prob,rprob)#randomize ties
        y=y[op]
        w=w[op]
        cw=cumsum(w)
        w1=w[y==1]
        cw1=cumsum(w1)
        wauc = log(sum(w1*(cw[y==1]-cw1)))
        sumw1 = cw1[length(cw1)]
        sumw2  = cw[length(cw)] - sumw1
        exp(wauc - log(sumw1) - log(sumw2))
    }
}

#' compute gradient for cox model
#'
#' Compute the gradient of the partial likelihood at a particular fit
#'
#' Compute a gradient vector at the fitted vector for the log partial likelihood.
#' This is like a residual vector, and useful for manual screening of predictors for \code{glmnet}
#' in applications where \code{p} is very large (as in GWAS). Uses the Breslow approach to ties
#'
#' @aliases coxgrad
#' @param f fit vector
#' @param time time vector (can have ties)
#' @param d death/censoring indicator 1/0
#' @param w observation weights (default equal)
#' @param eps (default 0.00001) Breaks ties between death and censoring by making death times \code{eps} earlier
#' @return a single gradient vector the same length as \code{f}
#' @author Trevor Hastie\cr Maintainer: Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{coxnet.deviance}
#' @keywords Cox model
#'
#' @export coxgrad

coxgrad=function(f,time,d,w,eps=0.00001){
### f is fitted function from glmnet at a particular lambda
### time is death or censoring time
### d is death indicator; d=0 means censored, d=1 means death
### w is a weight vector of non-negative weights, which will be normalized to sum to 1
    if(missing(w))w=rep(1,length(f))
    w=w/sum(w)
    f=scale(f,TRUE,FALSE)#center f so exponents are not too large
    time=time-d*eps#break ties between death times and non death times, leaving tied death times tied
    o=order(time)
    ef=exp(f)[o]
    time=time[o]
    d=d[o]
    w=w[o]
    rskden=rev(cumsum(rev(ef*w))) ##reverse order inside;last guy is in all the risk sets
### See if there are dups in death times
    dups=fid(time[d==1],seq(length(d))[d==1])
    dd=d
    ww=w
### next code replaces each sequence of tied death indicators by a new
### sequence where only the first is a 1 and the rest are zero. This
### makes the accounting in the following step work properly we also
### sums the weights in each of the tied death sets, and assign that
### weight to the first
    if(!is.null(ties<-dups$index_ties)){
        dd[unlist(ties)]=0
        dd[dups$index_first]=1
        wsum=sapply(ties,function(i,w)sum(w[i]),ww)
        tie1=sapply(ties,function(i)i[1])
        ww[tie1]=wsum
        }
### Get counts over risk sets at each death time
    rskcount=cumsum(dd)#this says how many of the risk sets each observation is in; 0 is none
### We now form partial sums of the 1/den just at the risk sets
    rskdeninv=cumsum((ww/rskden)[dd==1])
### pad with a zero, so we can index it
    rskdeninv=c(0,rskdeninv)
### compute gradient for each obs
    grad=(d-rskdeninv[rskcount+1]*ef)*w
    grad[o]=grad
    grad
    }

fid=function(x,index){
### Input:
### x is a sorted vector of death times
### index is vector of indices of this set
### Output:
### index of first member of every death set as they appear in sorted list
### list of ties for each element of index, in the case of two or more ties;
    ## if no ties, this list is NULL
    idup=duplicated(x)
    if(!any(idup)) list(index_first=index,index_ties=NULL)
    else{
        ndup=!idup
        xu=x[ndup]# first death times
        index_first=index[ndup]
        ities=match(x,xu)
        index_ties=split(index,ities)
        nties=sapply(index_ties,length)
        list(index_first=index_first,index_ties=index_ties[nties>1])
        }
    }


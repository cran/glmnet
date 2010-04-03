response.coxnet=function(y){
    if(!is.matrix(y)||!all(match(c("time","status"),dimnames(y)[[2]],0)))stop("Cox model requires a matrix with columns 'time' (>0) and 'status'  (binary) as a response; a 'Surv' object suffices",call.=FALSE)
  ty=as.double(y[,"time"])
  tevent=as.double(y[,"status"])
  if(any(ty<=0))stop("negative event times encountered;  not permitted for Cox family")
list(time=ty,event=tevent)
  }

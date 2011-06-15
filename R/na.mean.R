na.mean=function(x){
  ### replace missing entries with mean for that column.
  namean=function(x){
    nas=is.na(x)
    if(any(nas))  x[nas]=mean(x,na.rm=TRUE)
    x
  }
  apply(x,2,namean)
}

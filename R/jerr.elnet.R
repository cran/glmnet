jerr.elnet=function(n){
  if(n<7777)msg="Memory allocation error; contact package maintainer"
  else if(n==7777)msg="All used predictors have zero variance"
  else if(n==10000)msg="All penalty factors are <= 0"
  else msg="Unknown error"
  list(n,fatal=TRUE,msg)
}

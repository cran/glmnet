print.glmnet=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")

  print(cbind(Df=x$df,"%Dev"=signif(x$dev.ratio,digits),Lambda=signif(x$lambda,digits)))
    }

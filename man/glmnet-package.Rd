\name{glmnet-package}
\alias{glmnet-package}
\docType{package}
\title{
Elastic net model paths for some generalized linear models
}
\description{
  This package fits lasso and elastic-net model paths for regression, logistic
  and multinomial regression using coordinate descent. The algorithm is
  extremely fast, and exploits sparsity in the input x matrix where it
  exists. A variety of predictions can be made from the fitted models.
}
\details{
\tabular{ll}{
Package: \tab glmnet\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-05-14\cr
License: \tab What license is it under?\cr
}
Very simple to use. Accepts \code{x,y} data for regression models, and
produces the regularization path over a grid of values for the tuning
parameter \code{lambda}. Only 5 functions:
\code{glmnet}\cr
\code{predict.glmnet}\cr
\code{plot.glmnet}\cr
\code{print.glmnet}\cr
\code{coef.glmnet}
}
\author{Jerome Friedman, Trevor Hastie and Rob Tibshirani\cr
Maintainer: Trevor Hastie  <hastie@stanford.edu>
}
\references{Friedman, J., Hastie, T. and Tibshirani, R. (2008)
  \emph{Regularization Paths for Generalized Linear Models via Coordinate
    Descent},   \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}\cr
  \emph{Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr
  Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011)
  \emph{Regularization Paths for Cox's Proportional Hazards Model via
    Coordinate Descent, Journal of Statistical Software, Vol. 39(5)
    1-13}\cr
  \url{http://www.jstatsoft.org/v39/i05/}
}
\keyword{models}
\keyword{regression}
\keyword{ package }
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
g2=sample(1:2,100,replace=TRUE)
g4=sample(1:4,100,replace=TRUE)
fit1=glmnet(x,y)
predict(fit1,newx=x[1:5,],s=c(0.01,0.005))
predict(fit1,type="coef")
plot(fit1,xvar="lambda")
fit2=glmnet(x,g2,family="binomial")
predict(fit2,type="response",newx=x[2:5,])
predict(fit2,type="nonzero")
fit3=glmnet(x,g4,family="multinomial")
predict(fit3,newx=x[1:3,],type="response",s=0.01)
}

\name{beta_CVX}
\alias{x}
\alias{y}
\alias{beta_CVX}
\title{Simulated data for the glmnet vignette}
\description{
Simple simulated data, used to demonstrate the features of glmnet
}
\usage{
data(BinomialExample)
data(CVXResults)
data(CoxExample)
data(MultiGaussianExample)
data(MultinomialExample)
data(PoissonExample)
data(QuickStartExample)
data(SparseExample)
}
\format{
  Data objects used to demonstrate features in the glmnet vignette
}
\details{
Theses datasets are artificial, and ere used to test out some of the features
of glmnet.
}
\examples{
data(QuickStartExample)
glmnet(x,y)
}
\keyword{datasets}

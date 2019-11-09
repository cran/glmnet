#' Generate multinomial samples from a probability matrix
#'
#' Generate multinomial samples
#'
#' Simple function that calls the \code{rmultinom} function. It generates a class label
#' for each row of its input matrix of class probabilities.
#'
#' @param p matrix of probabilities, with number of columns the number of classes
#' @return a vector of class memberships
#' @author Trevor Hastie \cr Maintainer: Trevor Hastie <hastie@@stanford.edu>
#' @export rmult
rmult <-
  function(p){
    x=t(apply(p,1,function(x)rmultinom(1,1,x)))
    x=x%*%seq(ncol(p))
    drop(x)
  }


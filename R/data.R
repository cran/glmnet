#' Synthetic dataset with binary response
#'
#' Randomly generated data for binomial regression example.
#'
#' @name BinomialExample
#' @docType data
#' @usage data(BinomialExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{x}{100 by 30 matrix of numeric values.}
#'   \item{y}{Numeric vector of length 100 containing 44 zeros and 56 ones.}
#' }
"BinomialExample"

#' Synthetic dataset with right-censored survival response
#'
#' Randomly generated data for Cox regression example.
#'
#' @name CoxExample
#' @docType data
#' @usage data(CoxExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{x}{1,000 by 30 matrix of numeric values.}
#'   \item{y}{1,000 by 2 matrix with column names "time" and "status". The
#'     first column consists of positive numbers representing time to event,
#'     while the second column represents the status indicator
#'     (0=right-censored, 1=observed).}
#' }
"CoxExample"

#' Synthetic dataset with multiple Gaussian responses
#'
#' Randomly generated data for multi-response Gaussian regression example.
#'
#' @name MultiGaussianExample
#' @docType data
#' @usage data(MultiGaussianExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{x}{100 by 20 matrix of numeric values.}
#'   \item{y}{100 by 4 matrix of numeric values, each column representing
#'     one response vector.}
#' }
"MultiGaussianExample"

#' Synthetic dataset with multinomial response
#'
#' Randomly generated data for multinomial regression example.
#'
#' @name MultinomialExample
#' @docType data
#' @usage data(MultinomialExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{x}{500 by 30 matrix of numeric values.}
#'   \item{y}{Numeric vector of length 500 containing 142 ones, 174 twos
#'     and 184 threes.}
#' }
"MultinomialExample"

#' Synthetic dataset with count response
#'
#' Randomly generated data for Poisson regression example.
#'
#' @name PoissonExample
#' @docType data
#' @usage data(PoissonExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{x}{500 by 20 matrix of numeric values.}
#'   \item{y}{Numeric vector of length 500 consisting of non-negative
#'     integers.}
#' }
"PoissonExample"

#' Synthetic dataset with Gaussian response
#'
#' Randomly generated data for Gaussian regression example.
#'
#' @name QuickStartExample
#' @docType data
#' @usage data(QuickStartExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{x}{100 by 20 matrix of numeric values.}
#'   \item{y}{Numeric vector of length 100.}
#' }
"QuickStartExample"

#' Synthetic dataset with sparse design matrix
#'
#' Randomly generated data for Gaussian regression example with the
#' design matrix x being in sparse matrix format.
#'
#' @name SparseExample
#' @docType data
#' @usage data(SparseExample)
#' @keywords data
#' @format List containing the following elements:
#' \describe{
#'   \item{x}{100 by 20 matrix of numeric values. x is in sparse matrix
#'     format, having class "dgCMatrix".}
#'   \item{y}{Numeric vector of length 100.}
#' }
"SparseExample"

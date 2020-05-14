#' Display the names of the measures used in CV for different "glmnet" families
#'
#' Produces a list of names of measures
#'
#' Try it and see. A very simple function to provide information
#'
#' @param family If a "glmnet" family is supplied, a list of the names of
#' measures available for that family are produced. Default is "all", in which
#' case the names of measures for all families are produced.
#' @author Trevor Hastie\cr Maintainer: Trevor Hastie
#' \email{hastie@@stanford.edu}
#' @seealso \code{cv.glmnet} and \code{assess.glmnet}.
#' @keywords models
#' @export glmnet.measures
glmnet.measures=function(family=c("all","gaussian", "binomial", "poisson",
                                  "multinomial", "cox", "mgaussian","GLM")){
    family=match.arg(family)
    type.measures = c("mse", "deviance", "class", "auc", "mae",
                      "C")
    family.ch = list( gaussian = c(1, 5), binomial = c(2,
                                                          3, 4, 1, 5), poisson = c(2, 1, 5), cox = c(2, 6),
                     multinomial = c(2, 3, 1, 5), mgaussian = c(1, 5),GLM=c(2,1,5))
    if(family=="all")lapply(family.ch,function(x)type.measures[x])
    else         type.measures[family.ch[[family]]]

        }

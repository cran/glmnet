#' assess performance of a 'glmnet' object using test data.
#'
#' Given a test set, produce summary performance measures for the glmnet
#' model(s)
#'
#' \code{assess.glmnet} produces all the different performance measures
#' provided by \code{cv.glmnet} for each of the families. A single vector, or a
#' matrix of predictions can be provided, or fitted model objects or CV
#' objects. In the case when the predictions are still to be made, the
#' \code{...} arguments allow, for example, 'offsets' and other prediction
#' parameters such as values for 'gamma' for 'relaxed' fits.  \code{roc.glmnet}
#' produces for a single vector a two column matrix with columns TPR and FPR
#' (true positive rate and false positive rate). This object can be plotted to
#' produce an ROC curve. If more than one predictions are called for, then a
#' list of such matrices is produced.  \code{confusion.glmnet} produces a
#' confusion matrix tabulating the classification results. Again, a single
#' table or a list, with a print method.
#'
#' @aliases assess.glmnet confusion.glmnet roc.glmnet
#' @param object Fitted \code{"glmnet"} or \code{"cv.glmnet"}, \code{"relaxed"}
#' or \code{"cv.relaxed"} object, OR a matrix of predictions (for
#' \code{roc.glmnet} or \code{assess.glmnet}). For \code{roc.glmnet} the model
#' must be a 'binomial', and for \code{confusion.glmnet} must be either
#' 'binomial' or 'multinomial'
#' @param newx If predictions are to made, these are the 'x' values. Required
#' for \code{confusion.glmnet}
#' @param newy required argument for all functions; the new response values
#' @param weights For observation weights for the test observations
#' @param family The family of the model, in case predictions are passed in as
#' 'object'
#' @param ... additional arguments to \code{predict.glmnet} when "object" is a
#' "glmnet" fit, and predictions must be made to produce the statistics.
#' @return \code{assess.glmnet} produces a list of vectors of measures.
#' \code{roc.glmnet} a list of 'roc' two-column matrices, and
#' \code{confusion.glmnet} a list of tables. If a single prediction is
#' provided, or predictions are made from a CV object, the latter two drop the
#' list status and produce a single matrix or table.
#' @author Trevor Hastie and Rob Tibshirani\cr Maintainer: Trevor Hastie
#' <hastie@@stanford.edu>
#' @seealso \code{cv.glmnet}, \code{glmnet.measures} and \code{vignette("relax",package="glmnet")}
#' @keywords models classification
#' @examples
#'
#' data(QuickStartExample)
#' x <- QuickStartExample$x; y <- QuickStartExample$y
#' set.seed(11)
#' train = sample(seq(length(y)),70,replace=FALSE)
#' fit1 = glmnet(x[train,], y[train])
#' assess.glmnet(fit1, newx = x[-train,], newy = y[-train])
#' preds = predict(fit1, newx = x[-train, ], s = c(1, 0.25))
#' assess.glmnet(preds, newy = y[-train], family = "gaussian")
#' fit1c = cv.glmnet(x, y, keep = TRUE)
#' fit1a = assess.glmnet(fit1c$fit.preval, newy=y,family="gaussian")
#' plot(fit1c$lambda, log="x",fit1a$mae,xlab="Log Lambda",ylab="Mean Absolute Error")
#' abline(v=fit1c$lambda.min, lty=2, col="red")
#' data(BinomialExample)
#' x <- BinomialExample$x; y <- BinomialExample$y
#' fit2 = glmnet(x[train,], y[train], family = "binomial")
#' assess.glmnet(fit2,newx = x[-train,], newy=y[-train], s=0.1)
#' plot(roc.glmnet(fit2, newx = x[-train,], newy=y[-train])[[10]])
#' fit2c = cv.glmnet(x, y, family = "binomial", keep=TRUE)
#' idmin = match(fit2c$lambda.min, fit2c$lambda)
#' plot(roc.glmnet(fit2c$fit.preval, newy = y)[[idmin]])
#' data(MultinomialExample)
#' x <- MultinomialExample$x; y <- MultinomialExample$y
#' set.seed(103)
#' train = sample(seq(length(y)),100,replace=FALSE)
#' fit3 = glmnet(x[train,], y[train], family = "multinomial")
#' confusion.glmnet(fit3, newx = x[-train, ], newy = y[-train], s = 0.01)
#' fit3c = cv.glmnet(x, y, family = "multinomial", type.measure="class", keep=TRUE)
#' idmin = match(fit3c$lambda.min, fit3c$lambda)
#' confusion.glmnet(fit3c$fit.preval, newy = y, family="multinomial")[[idmin]]
#'
#' @export assess.glmnet
assess.glmnet=function(object,newx=NULL,newy,weights=NULL,
                       family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),
                       ...){
### object must be either a glmnet or cv.glmnet object, or else a matrix/vector of predictions of a glmnet model  (not on the mean scale, but natural parameter scale)
### It must either be a matrix or an array. The last dimension can be 1.
    oclass=grep("glmnet",class(object))
    if(length(oclass)){# either a glmnet or cv.glmnet
        predmat=predict(object,newx=newx,...)
        fam=family(object)
    }
    else{
        if(is.character(family))fam=match.arg(family)
        else {
            if(is.function(family))family=family()
            if(!inherits(family,"family"))
                stop("Invalid family argument; must be either character, function or family object")
            fam=family
            }

        predmat=object
        if(is.vector(predmat))predmat=as.matrix(predmat)
    }
    famtype=fam
    attr(predmat,"family")=fam
    if(inherits(fam,"family"))famtype="GLM"
    type.measures = glmnet.measures(famtype)
    y = drop(newy)
    dimy=dim(y)
    nrowy=ifelse(is.null(dimy),length(y),dimy[1])
    if (is.null(weights))
        weights = rep(1, nrowy)
    subclass=switch(famtype,gaussian="elnet",binomial="lognet",poisson="fishnet",multinomial="multnet",
                    cox="coxnet",mgaussian="mrelnet",GLM="glmnetfit")
    fun = paste("cv", subclass, sep = ".")
    if(subclass=="coxnet")fun="assess.coxnet"
    outlist=as.list(type.measures)
    names(outlist)=type.measures
    for(type.measure in type.measures){
    type.measure=cvtype(type.measure,subclass)
    teststuff = do.call(fun, list(predmat,y,type.measure,weights,foldid=rep(1,nrowy),grouped=FALSE))
    out=drop(with(teststuff, apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)))
    attr(out,"measure")=names(type.measure)
    outlist[[type.measure]]=out
    }
    outlist
}


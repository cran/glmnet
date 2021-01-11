#' convert a data frame to a data matrix with one-hot encoding
#'
#' Converts a data frame to a data matrix suitable for input to \code{glmnet}.
#' Factors are converted to dummy matrices via "one-hot" encoding. Options deal
#' with missing values and sparsity.
#'
#' The main function is to convert factors to dummy matrices via "one-hot"
#' encoding. Having the 'train' and 'test' data present is useful if some
#' factor levels are missing in either. Since a factor with k levels leads to a
#' submatrix with 1/k entries zero, with large k the \code{sparse=TRUE} option
#' can be helpful; a large matrix will be returned, but stored in sparse matrix
#' format. Finally, the function can deal with missing data. The current
#' version has the option to replace missing observations with the mean from
#' the training data. For dummy submatrices, these are the mean proportions at
#' each level.
#'
#' @param train Required argument. A dataframe consisting of vectors, matrices
#' and factors
#' @param test Optional argument. A dataframe matching 'train' for use as
#' testing data
#' @param na.impute Logical, default \code{FALSE}. If \code{TRUE}, missing
#' values for any column in the resultant 'x' matrix are replaced by the means
#' of the nonmissing values derived from 'train'
#' @param sparse Logical, default \code{FALSE}. If \code{TRUE} then the
#' returned matrice(s) are converted to matrices of class "CsparseMatrix".
#' Useful if some factors have a large number of levels, resulting in very big
#' matrices, mostly zero
#' @param ... additional arguments, currently unused
#' @return If only 'train' was provided, the function returns a matrix 'x'. If
#' missing values were imputed, this matrix has an attribute containing its
#' column means (before imputation). If 'test' was provided as well, a list
#' with two components is returned: 'x' and 'xtest'.
#' @author Trevor Hastie\cr Maintainer: Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{glmnet}
#' @keywords models
#' @examples
#'
#' set.seed(101)
#' ### Single data frame
#' X = matrix(rnorm(20), 10, 2)
#' X3 = sample(letters[1:3], 10, replace = TRUE)
#' X4 = sample(LETTERS[1:3], 10, replace = TRUE)
#' df = data.frame(X, X3, X4)
#' makeX(df)
#' makeX(df, sparse = TRUE)
#'
#' ### Single data freame with missing values
#' Xn = X
#' Xn[3, 1] = NA
#' Xn[5, 2] = NA
#' X3n = X3
#' X3n[6] = NA
#' X4n = X4
#' X4n[9] = NA
#' dfn = data.frame(Xn, X3n, X4n)
#'
#' makeX(dfn)
#' makeX(dfn, sparse = TRUE)
#' makeX(dfn, na.impute = TRUE)
#' makeX(dfn, na.impute = TRUE, sparse = TRUE)
#'
#' ### Test data as well
#' X = matrix(rnorm(10), 5, 2)
#' X3 = sample(letters[1:3], 5, replace = TRUE)
#' X4 = sample(LETTERS[1:3], 5, replace = TRUE)
#' dft = data.frame(X, X3, X4)
#'
#' makeX(df, dft)
#' makeX(df, dft, sparse = TRUE)
#'
#' ### Missing data in test as well
#' Xn = X
#' Xn[3, 1] = NA
#' Xn[5, 2] = NA
#' X3n = X3
#' X3n[1] = NA
#' X4n = X4
#' X4n[2] = NA
#' dftn = data.frame(Xn, X3n, X4n)
#'
#' makeX(dfn, dftn)
#' makeX(dfn, dftn, sparse = TRUE)
#' makeX(dfn, dftn, na.impute = TRUE)
#' makeX(dfn, dftn, sparse = TRUE, na.impute = TRUE)
#'
#' @export makeX
makeX=function(train,test=NULL, na.impute=FALSE,sparse=FALSE,...){
    df=train
    istest=!is.null(test)
    ntr=nrow(train)
    if(istest){
        nte=nrow(test)
        df=rbind(df,test)
    }
    ### bug fix because of change of default behavior of data.frame
    ## check if any character columns
    classes = sapply(df,class)
    if (any(classes == "character"))
        df <- as.data.frame(unclass(df), stringsAsFactors = TRUE)
    ###
    x=prepareX(df,sparse=sparse)
    if(na.impute){
        xbar=colMeans(x[seq(ntr),],na.rm=TRUE)
        x=na.replace(x,xbar)
    }
    if(istest){
        xt=x[seq(ntr+1,ntr+nte),]
        x=x[seq(ntr),]
        }
        if(na.impute)attr(x,"means")=xbar
    if(istest)list(x=x,xtest=xt) else x
}

#' @importFrom stats contrasts model.matrix na.pass
#' @importFrom Matrix sparse.model.matrix
#' @export prepareX
prepareX <-
    function(df,sparse=FALSE,...){
        if(!inherits(df,"data.frame"))stop("first argument must be of class `data.frame`")
        whichfac=sapply(df,inherits,"factor")
### bug fix if factor has one level only
        ## Cannot get contrasts to do the job here, so do it manually
        df_level <- lapply(df[,whichfac,drop = FALSE], levels)
        nlevels <- sapply(df_level,length)
        which <- nlevels == 1
        if(any(which)){
            whichfac1 = whichfac
            whichfac1[whichfac] = which
            dn <- names(df)
            cx <- as.matrix(df[,whichfac1],drop=FALSE)
            unimat <- array(1,dim=dim(cx))
            cnames <- paste0(dn[whichfac1],unlist(df_level[which]))
            unimat[is.na(cx)] <- NA
            df[,whichfac1] <- data.frame(unimat)
            dn[whichfac1] <- cnames
            names(df) <- dn
            whichfac[whichfac] = !which
            warning(call. = FALSE,paste("Column(s) ",paste(cnames,collapse=", "), "are all 1, due to factors with a single level"))
        }
###
        oldna=options()$na.action
        cna=as.character(substitute(na.action))
        options(na.action=na.pass)
        on.exit(options(na.action=oldna))
        if(any(whichfac))
            ctr=lapply(df[,whichfac,drop=FALSE], contrasts,contrast=FALSE)
        else ctr=NULL
        if(sparse){
            m=sparse.model.matrix(~.-1,data=df,contrasts.arg=ctr,...)
            m=na_sparse_fix(m,names(df)) # sparse.model.matrix is faulty
            }
           else m = model.matrix(~.-1,data=df,contrasts.arg=ctr,...)
        if(any(whichfac))attr(m,"contrasts")=NULL
        attr(m,"assign")=NULL
        m
    }

#' @export na_sparse_fix
na_sparse_fix=function(x,dfnames){
    a=attributes(x)
    ac=a$contrasts
    as=a$assign
    if(is.null(ac))return(x)
    acn=names(ac)
    whichn=match(acn,dfnames)
    for(i in whichn){
        xi=x[,as==i]
        rowtot=rowSums(xi)
        if(sum(rowtot)<length(rowtot)){# Nas present
            x[rowtot==0,as==i]=NA
        }
    }
    x
}


#' Replace the missing entries in a matrix columnwise with the entries in a
#' supplied vector
#'
#' Missing entries in any given column of the matrix are replaced by the column
#' means or the values in a supplied vector.
#'
#' This is a simple imputation scheme. This function is called by \code{makeX}
#' if the \code{na.impute=TRUE} option is used, but of course can be used on
#' its own. If 'x' is sparse, the result is sparse, and the replacements are
#' done so as to maintain sparsity.
#'
#' @param x A matrix with potentially missing values, and also potentially in
#' sparse matrix format (i.e. inherits from "sparseMatrix")
#' @param m Optional argument. A vector of values used to replace the missing
#' entries, columnwise. If missing, the column means of 'x' are used
#' @return A version of 'x' is returned with the missing values replaced.
#' @author Trevor Hastie\cr Maintainer: Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{makeX} and \code{glmnet}
#' @keywords models
#' @examples
#'
#' set.seed(101)
#' ### Single data frame
#' X = matrix(rnorm(20), 10, 2)
#' X[3, 1] = NA
#' X[5, 2] = NA
#' X3 = sample(letters[1:3], 10, replace = TRUE)
#' X3[6] = NA
#' X4 = sample(LETTERS[1:3], 10, replace = TRUE)
#' X4[9] = NA
#' dfn = data.frame(X, X3, X4)
#'
#' x = makeX(dfn)
#' m = rowSums(x, na.rm = TRUE)
#' na.replace(x, m)
#'
#' x = makeX(dfn, sparse = TRUE)
#' na.replace(x, m)
#'
#' @export na.replace
na.replace=function(x,m=rowSums(x,na.rm=TRUE)){
    if(inherits(x,"sparseMatrix")){
        x=as(x,"CsparseMatrix")
        ccount=diff(x@p)
        cindex=rep(1:length(ccount),ccount)
        nas=is.na(x@x)
        if(any(nas))x@x[nas]=m[cindex[nas]]
        }
    else{
            d=dim(x)
            cols=rep(1:d[2],rep(d[1],d[2]))
            nas=is.na(x)
            if(any(nas))x[nas]=m[cols[nas]]
            }
    x
}

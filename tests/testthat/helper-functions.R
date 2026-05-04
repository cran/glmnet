# Helper functions for glmnet tests
# This file is auto-sourced by testthat before running tests

library(survival)
library(nnet)

## .control_with(): build a fully resolved control list for direct calls
## to glmnet.path / glmnet.fit / elnet.fit, which no longer accept scalar
## thresh / maxit / trace.it args in their signatures. Merges user
## overrides onto the current session defaults from glmnet.control(),
## validating key names against the canonical .VALID_CONTROL_KEYS.
##
## Usage in tests:
##   glmnet.path(x, y, ..., control = .control_with(thresh = 1e-8))
##   glmnet.fit (x, y, ..., control = .control_with(thresh = 1e-10, maxit = 500))
##
## Test-only helper; not exported by the package.
.control_with <- function(...) {
    args <- list(...)
    if (length(args) && (is.null(names(args)) || any(names(args) == "")))
        stop(".control_with() requires named arguments", call. = FALSE)
    valid <- glmnet:::.VALID_CONTROL_KEYS
    bad <- setdiff(names(args), valid)
    if (length(bad))
        stop("invalid control key(s) in .control_with(): ",
             paste(bad, collapse = ", "),
             "\n  valid keys: ", paste(valid, collapse = ", "),
             call. = FALSE)
    modifyList(glmnet.control(), args)
}

# Remove $call from a glmnet object (and nested objects like $glmnet.fit)
# so that expect_equal comparisons are not affected by changes in how
# arguments are passed (e.g., thresh= vs control=list(thresh=))
drop_call <- function(obj) {
  obj$call <- NULL
  if (!is.null(obj$glmnet.fit)) {
    obj$glmnet.fit$call <- NULL
  }
  if (!is.null(obj$relaxed)) {
    obj$relaxed$call <- NULL
  }
  obj
}

# Test mode control: "test" (default), "save", or "overwrite"
# - "test": run tests against saved results in saved_results/
# - "save": save results to timestamped folder saved_results/YYYY-MM-DD:HH:MM:SS/
# - "overwrite": overwrite saved results in saved_results/
# - any other value: skip everything
test_modes <- c("test", "save", "overwrite")
test_mode <- "test"

# Check if test_mode is valid - if not, tests should skip everything
# Can pass a different test_mode to override the default
valid_test_mode <- function(mode = test_mode) {
    if (mode %in% test_modes) {
        return(TRUE)
    }
    message(sprintf("Skipping tests: test_mode '%s' is not valid. Expected one of: %s",
                    mode, paste(test_modes, collapse = ", ")))
    FALSE
}

# Save test results based on test_mode
# - "test": does nothing (returns invisibly)
# - "save": saves to timestamped folder saved_results/YYYY-MM-DD:HH:MM:SS/
# - "overwrite": overwrites saved results in saved_results/
save_test_results <- function(fit.list, filename, mode = test_mode) {
    if (mode == "test") {
        return(invisible(NULL))
    }
    if (mode == "save") {
        timestamp_dir <- format(Sys.time(), "%Y-%m-%d:%H:%M:%S")
        save_dir <- file.path("saved_results", timestamp_dir)
        if (!dir.exists(save_dir)) {
            dir.create(save_dir, recursive = TRUE)
        }
        save_path <- file.path(save_dir, filename)
    } else if (mode == "overwrite") {
        save_path <- file.path("saved_results", filename)
    } else {
        return(invisible(NULL))
    }
    saveRDS(fit.list, file = save_path)
    message("Saved results to: ", save_path)
    invisible(save_path)
}

## Utility functions from preamble.R
expectations <- function(file) {
    expected <- readRDS(file)
    for (x in names(expected)) {
        cat(sprintf("expect_identical(%s, expected$%s)\n", x, x))
    }
    invisible()
}

ck <- function(x, y, item) {
    nx <- names(x)
    ny <- names(y)
    if (identical(sort(names(x)), sort(c(names(y), item)))) {
        x[[item]] <- NULL
        expect_equal(x, y)
    } else {
        FALSE
    }
}

runtest <- function(fname = "elnet", saved = TRUE, rootdir = ".") {
    if (saved) cat("Testing", fname, "...\n") else cat("Saving", fname, "...\n")
    sfname <- paste(rootdir, "/test_", fname, ".R", sep = "")
    saved <<- saved
    source(sfname)
    sobfile <- paste("Benchmark/save.", fname, ".RData", sep = "")
    if (!saved) {
        assign("save.list", fit.list)
        save(save.list, file = sobfile)
    } else {
        load(sobfile)
        snames <- names(fit.list)
        sanames <- names(save.list)
        for (i in snames) {
            if (!match(i, sanames, FALSE)) cat(i, " is new; not in saved objects\n")
            else {
                j1 <- all.equal(fit.list[[i]], save.list[[i]])
                cat(i, isTRUE(j1), "\n")
                if (!isTRUE(j1)) print(j1)
            }
        }
        invisible(fit.list)
    }
}

enlist <- function (...) {
    result <- list(...)
    if ((nargs() == 1) & is.character(n <- result[[1]])) {
        result <- as.list(seq(n))
        names(result) <- n
        for (i in n) result[[i]] <- get(i)
    }
    else {
        n <- sys.call()
        n <- as.character(n)[-1]
        if (!is.null(n2 <- names(result))) {
            which <- n2 != ""
            n[which] <- n2[which]
        }
        names(result) <- n
    }
    result
}

MakeYMat  <-  function(y) {
    ymat <- matrix(0, nrow=length(y), ncol=length(unique(y)))
    for(i in 1:ncol(ymat)) ymat[y==i,i] <- 1
    return(ymat)
}

norm2  <- function(x, ...) {
    Matrix::norm(x, type = "2", ...)
}

kkt.multnet <- function(x,y,standardize=TRUE,intercept=TRUE,which.lambda=15,offset=FALSE,digits=7,type.multinomial="ungrouped",...) {

    family="multinomial"
    kpass=function(mat)(sum(apply(mat,1,prod))==0)&all(mat[,2]<=0)
    if(missing(offset))
        fit=glmnet(x,y,intercept=intercept,standardize=standardize,family=family,type.multinomial=type.multinomial,...)
    else
        fit=glmnet(x,y,intercept=intercept,standardize=standardize,family=family,type.multinomial=type.multinomial,offset=offset,...)
    n = nrow(x)
    sdj=function(x)drop(sqrt(var(x)*(length(x)-1)/length(x)))
    if(standardize)sdx=apply(x,2,sdj)else sdx=rep(1,ncol(x))
    xsd=scale(x,intercept,sdx)
    lambda0=fit$lambda[which.lambda]
    ## print(coef(fit,s=lambda0))
    if(missing(offset))  pred=predict(fit,x,s=lambda0,type="response")
    else  pred=predict(fit,x,s=lambda0,type="response",newoffset=offset)
    resid=y-drop(pred)
    cmat=sapply(coef(fit,s=lambda0),as.matrix)
    if(type.multinomial=="ungrouped"){
        kmat=round(abs(t(xsd)%*%resid/n)-lambda0,digits)
        mat=array(rbind(cmat,rbind(-!intercept,kmat)),c(nrow(cmat),2,ncol(cmat)))
        kkt.pass= apply(mat,3,kpass)
    }
    else{
        nrm=apply(cmat,1,norm2)
        gmat=t(xsd)%*%resid/n
        gnorm=apply(gmat,1,norm2)
        kvec=round(gnorm-lambda0,digits)
        mat=cbind(nrm,c(-!intercept,kvec))
        kkt.pass=kpass(mat)
    }
    attr(mat,"kkt.pass")=kkt.pass
    mat
}


kkt.glmnet <- function(x,y,standardize=TRUE,intercept=TRUE,which.lambda=15,family="gaussian",offset=FALSE,digits=7,...) {
    if(missing(offset))
        fit=glmnet(x,y,intercept=intercept,standardize=standardize,family=family,...)
    else
        fit=glmnet(x,y,intercept=intercept,standardize=standardize,family=family,offset=offset,...)
    n = nrow(x)
    sdj=function(x)drop(sqrt(var(x)*(length(x)-1)/length(x)))
    if(standardize)sdx=apply(x,2,sdj)else sdx=rep(1,ncol(x))
    xsd=scale(x,intercept,sdx)
    if(family=="gaussian"){sdy=sdj(y);ysd=y/sdy}  else {sdy=1;ysd=y}
    lambda0=fit$lambda[which.lambda]
    ## print(coef(fit,s=lambda0))
    if(missing(offset))  pred=predict(fit,x,s=lambda0,type="response")/sdy
    else  pred=predict(fit,x,s=lambda0,type="response",newoffset=offset)/sdy

    resid=ysd-pred
    mat=cbind(as.matrix(coef(fit,s=lambda0)),rbind(-!intercept,round(abs(t(xsd)%*%resid/n)-lambda0/sdy,digits)))
    kkt.pass= (sum(apply(mat,1,prod))==0)&all(mat[,2]<=0)
    attr(kkt.pass,"kkt.mat")=mat
    kkt.pass
}

sysinfo  <- Sys.info()

# Helper function for approximate equality using relative and absolute tolerance
# For scalar comparison
approx_equal <- function(a, b, rtol = 1e-10, atol = 0.0) {
    abs(a - b) <= atol + rtol * max(abs(a), abs(b))
}

# For vector comparison - returns TRUE if all elements are approximately equal
all_approx_equal <- function(a, b, rtol = 1e-10, atol = 0.0) {
    all(abs(a - b) <= atol + rtol * pmax(abs(a), abs(b)))
}

# Get max relative/absolute difference for diagnostics
max_diff <- function(a, b) {
    list(
        abs = max(abs(a - b)),
        rel = max(abs(a - b) / pmax(abs(a), abs(b), 1e-100))
    )
}

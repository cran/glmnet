## ----include=FALSE------------------------------------------------------------
# the code in this chunk enables us to truncate the print output for each
# chunk using the `out.lines` option
# save the built-in output hook
hook_output <- knitr::knit_hooks$get("output")

# set a new output hook to truncate text output
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
        
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("glmnet", repos = "https://cran.us.r-project.org")

## -----------------------------------------------------------------------------
library(glmnet)

## -----------------------------------------------------------------------------
data(QuickStartExample)

## -----------------------------------------------------------------------------
fit <- glmnet(x, y)

## -----------------------------------------------------------------------------
plot(fit)

## ----out.lines = 10-----------------------------------------------------------
print(fit)

## ----out.lines = 10-----------------------------------------------------------
coef(fit, s = 0.1)

## -----------------------------------------------------------------------------
set.seed(29)
nx <- matrix(rnorm(5 * 20), 5, 20)
predict(fit, newx = nx, s = c(0.1, 0.05))

## -----------------------------------------------------------------------------
cvfit <- cv.glmnet(x, y)

## -----------------------------------------------------------------------------
plot(cvfit)

## ----out.lines = 10-----------------------------------------------------------
cvfit$lambda.min
coef(cvfit, s = "lambda.min")

## -----------------------------------------------------------------------------
predict(cvfit, newx = x[1:5,], s = "lambda.min")

## -----------------------------------------------------------------------------
wts <-  c(rep(1,50), rep(2,50))
fit <- glmnet(x, y, alpha = 0.2, weights = wts, nlambda = 20)

## -----------------------------------------------------------------------------
print(fit)

## -----------------------------------------------------------------------------
fit <- glmnet(x, y)
any(fit$lambda == 0.5)  # 0.5 not in original lambda sequence
coef.apprx <- coef(fit, s = 0.5, exact = FALSE)
coef.exact <- coef(fit, s = 0.5, exact = TRUE, x=x, y=y)
cbind2(coef.exact[which(coef.exact != 0)], 
       coef.apprx[which(coef.apprx != 0)])

## -----------------------------------------------------------------------------
predict(fit, newx = x[1:5,], type = "response", s = 0.05)

## -----------------------------------------------------------------------------
plot(fit, xvar = "lambda", label = TRUE)

## -----------------------------------------------------------------------------
plot(fit, xvar = "dev", label = TRUE)

## -----------------------------------------------------------------------------
cvfit <- cv.glmnet(x, y, type.measure = "mse", nfolds = 20)

## -----------------------------------------------------------------------------
print(cvfit)

## ---- eval=FALSE--------------------------------------------------------------
#  library(doMC)
#  registerDoMC(cores = 2)
#  X <- matrix(rnorm(1e4 * 200), 1e4, 200)
#  Y <- rnorm(1e4)

## ---- eval=FALSE--------------------------------------------------------------
#  system.time(cv.glmnet(X, Y))

## ---- echo=FALSE--------------------------------------------------------------
structure(c(2.44, 0.08, 2.518, 0, 0), class = "proc_time", .Names = c("user.self",
"sys.self", "elapsed", "user.child", "sys.child"))

## ---- eval=FALSE--------------------------------------------------------------
#  system.time(cv.glmnet(X, Y, parallel = TRUE))

## ---- echo=FALSE--------------------------------------------------------------
structure(c(0.508999999999999, 0.057, 1.56699999999999, 1.941,
0.1), class = "proc_time", .Names = c("user.self", "sys.self",
"elapsed", "user.child", "sys.child"))

## ----out.lines = 10-----------------------------------------------------------
cvfit$lambda.min
predict(cvfit, newx = x[1:5,], s = "lambda.min")
coef(cvfit, s = "lambda.min")

## -----------------------------------------------------------------------------
foldid <- sample(1:10, size = length(y), replace = TRUE)
cv1  <- cv.glmnet(x, y, foldid = foldid, alpha = 1)
cv.5 <- cv.glmnet(x, y, foldid = foldid, alpha = 0.5)
cv0  <- cv.glmnet(x, y, foldid = foldid, alpha = 0)

## -----------------------------------------------------------------------------
par(mfrow = c(2,2))
plot(cv1); plot(cv.5); plot(cv0)
plot(log(cv1$lambda)   , cv1$cvm , pch = 19, col = "red",
     xlab = "log(Lambda)", ylab = cv1$name)
points(log(cv.5$lambda), cv.5$cvm, pch = 19, col = "grey")
points(log(cv0$lambda) , cv0$cvm , pch = 19, col = "blue")
legend("topleft", legend = c("alpha= 1", "alpha= .5", "alpha 0"),
       pch = 19, col = c("red","grey","blue"))

## -----------------------------------------------------------------------------
tfit <- glmnet(x, y, lower.limits = -0.7, upper.limits = 0.5)
plot(tfit)

## -----------------------------------------------------------------------------
p.fac <- rep(1, 20)
p.fac[c(1, 3, 5)] <- 0
pfit <- glmnet(x, y, penalty.factor = p.fac)
plot(pfit, label = TRUE)

## -----------------------------------------------------------------------------
data(MultiGaussianExample)
mfit <- glmnet(x, y, family = "mgaussian")

## -----------------------------------------------------------------------------
plot(mfit, xvar = "lambda", label = TRUE, type.coef = "2norm")

## -----------------------------------------------------------------------------
predict(mfit, newx = x[1:5,], s = c(0.1, 0.01))

## -----------------------------------------------------------------------------
data(BinomialExample)

## -----------------------------------------------------------------------------
fit <- glmnet(x, y, family = "binomial")

## -----------------------------------------------------------------------------
predict(fit, newx = x[1:5,], type = "class", s = c(0.05, 0.01))

## -----------------------------------------------------------------------------
cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "class")

## -----------------------------------------------------------------------------
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se

## -----------------------------------------------------------------------------
data(MultinomialExample)

## -----------------------------------------------------------------------------
fit <- glmnet(x, y, family = "multinomial", type.multinomial = "grouped")
plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")

## -----------------------------------------------------------------------------
cvfit <- cv.glmnet(x, y, family = "multinomial", type.multinomial = "grouped")
plot(cvfit)

## -----------------------------------------------------------------------------
predict(cvfit, newx = x[1:10,], s = "lambda.min", type = "class")

## -----------------------------------------------------------------------------
data(PoissonExample)

## -----------------------------------------------------------------------------
fit <- glmnet(x, y, family = "poisson")

## -----------------------------------------------------------------------------
plot(fit)

## ----out.lines = 7------------------------------------------------------------
coef(fit, s = 1)
predict(fit, newx = x[1:5,], type = "response", s = c(0.1,1))

## -----------------------------------------------------------------------------
cvfit <- cv.glmnet(x, y, family = "poisson")

## -----------------------------------------------------------------------------
data(BinomialExample)
itrain <- 1:70
fit <- glmnet(x[itrain, ], y[itrain], family = "binomial", nlambda = 5)
assess.glmnet(fit, newx = x[-itrain, ], newy = y[-itrain])

## ---- eval=FALSE--------------------------------------------------------------
#  pred <- predict(fit, newx = x[-itrain, ])
#  assess.glmnet(pred, newy = y[-itrain], family = "binomial")

## -----------------------------------------------------------------------------
glmnet.measures()

## ----out.lines = 11-----------------------------------------------------------
cfit <- cv.glmnet(x[itrain, ], y[itrain], family = "binomial", nlambda = 30)
assess.glmnet(cfit, newx = x[-itrain, ], newy = y[-itrain])

## ----out.lines = 11-----------------------------------------------------------
assess.glmnet(cfit, newx = x[-itrain, ],newy = y[-itrain], s = "lambda.min")

## ----out.lines = 11-----------------------------------------------------------
cfit <- cv.glmnet(x, y, family = "binomial", keep = TRUE, nlambda = 30)
assess.glmnet(cfit$fit.preval, newy = y, family = "binomial")

## -----------------------------------------------------------------------------
cfit <- cv.glmnet(x, y, family = "binomial", type.measure = "auc", 
                  keep = TRUE)
rocs <- roc.glmnet(cfit$fit.preval, newy = y)

## -----------------------------------------------------------------------------
best <- cvfit$index["min",]
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "red")

## -----------------------------------------------------------------------------
data(MultinomialExample)
set.seed(101)
itrain <- sample(1:500, 400, replace = FALSE)
cfit <- cv.glmnet(x[itrain, ], y[itrain], family = "multinomial")
cnf <- confusion.glmnet(cfit, newx = x[-itrain, ], newy = y[-itrain])

## -----------------------------------------------------------------------------
print(cnf)

## -----------------------------------------------------------------------------
cfit <- cv.glmnet(x, y, family = "multinomial", type = "class", keep = TRUE)
cnf <- confusion.glmnet(cfit$fit.preval, newy = y, family = "multinomial")
best <- cfit$index["min",]
print(cnf[[best]])

## -----------------------------------------------------------------------------
data(SparseExample)
class(x)

## -----------------------------------------------------------------------------
fit <- glmnet(x, y)

## -----------------------------------------------------------------------------
cvfit = cv.glmnet(x, y)
plot(cvfit)

## -----------------------------------------------------------------------------
i <- sample(1:5, size = 25, replace = TRUE)
j <- sample(1:20, size = 25, replace = TRUE)
x <- rnorm(25)
nx <- sparseMatrix(i = i, j = j, x = x, dims = c(5, 20))
predict(cvfit, newx = nx, s = "lambda.min")

## -----------------------------------------------------------------------------
data(BinomialExample)
fit <- bigGlm(x, y, family = "binomial", lower.limits = -1)
print(fit)

## -----------------------------------------------------------------------------
set.seed(101)
X <- matrix(rnorm(5), nrow = 5)
X2 <- sample(letters[1:3], 5, replace = TRUE)
X3 <- sample(LETTERS[1:3], 5, replace = TRUE)
df <- data.frame(X, X2, X3)
makeX(df)

## -----------------------------------------------------------------------------
makeX(df, sparse = TRUE)

## -----------------------------------------------------------------------------
Xn  <- X ; Xn[3,1] <- NA
X2n <- X2; X2n[1]  <- NA
X3n <- X3; X3n[5]  <- NA
dfn <- data.frame(Xn, X2n, X3n)
dfn
makeX(dfn)

## -----------------------------------------------------------------------------
makeX(dfn, na.impute = TRUE, sparse = TRUE)

## -----------------------------------------------------------------------------
set.seed(102)
X <- matrix(rnorm(5), nrow = 5)
X2 <- sample(letters[1:3], 5, replace = TRUE)
X3 <- sample(LETTERS[1:3], 5, replace = TRUE)
Xn  <- X ; Xn[5,1] <- NA
X2n <- X2; X2n[1]  <- NA
X3n <- X3; X3n[2]  <- NA
dftn <- data.frame(Xn, X2n, X3n)
dftn
makeX(dfn, dftn, sparse = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  fit <- glmnet(x, y, trace.it = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  fit <- cv.glmnet(x, y, trace.it = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  glmnet.control(itrace = 1)

## -----------------------------------------------------------------------------
data(QuickStartExample)
fit <- glmnet(x, y)
length(fit$lambda)  # number of lambda values fit

## -----------------------------------------------------------------------------
glmnet.control(fdev = 0.1)
fit <- glmnet(x, y)
length(fit$lambda)  # number of lambda values fit

## -----------------------------------------------------------------------------
glmnet.control(factory = TRUE)

## ----out.lines = 8------------------------------------------------------------
glmnet.control()

## ---- echo=FALSE--------------------------------------------------------------
data(QuickStartExample)

## -----------------------------------------------------------------------------
np <- dim(x); n <- np[1]; p <-np[2]

fit <- glmnet(x, y, intercept = F, standardize = F, 
              lambda = 8 / (2 * n), thresh = 1e-20)

## ----eval=FALSE---------------------------------------------------------------
#  beta_glmnet <- as.matrix(predict(fit, type = "coefficients")[-1,])

## -----------------------------------------------------------------------------
fit <- glmnet(x, y, intercept = F, standardize = F, thresh = 1e-20)
beta_glmnet <- as.matrix(predict(fit, s = 8 / (2 * n), 
                                 type = "coefficients", 
                                 exact = TRUE, x = x, y = y)[-1,])

## ---- eval=FALSE--------------------------------------------------------------
#  library(CVXR)
#  beta <- Variable(p)
#  loss <- sum((y-x%*%beta)^2)/(2*n)
#  lassoPenalty <- function(beta,lambda)lambda*p_norm(beta,1)
#  obj <- loss + lassoPenalty(beta, lambda = 8/(2*n))
#  prob <- Problem(Minimize(obj))
#  result <- solve(prob)
#  beta_CVX <- result$getValue(beta)

## -----------------------------------------------------------------------------
data(CVXResults)

## ---- message=FALSE-----------------------------------------------------------
library(lars)
fit_lars <- lars(x, y, type = "lasso", intercept = F, normalize = F)
beta_lars <- predict(fit_lars, s = 8 / 2, type = "coefficients", 
                     mode = "lambda")$coefficients

## -----------------------------------------------------------------------------
cmp <- round(cbind(beta_glmnet, beta_lars, beta_CVX), digits = 6)
colnames(cmp) <- c("beta_glmnet", "beta_lars", "beta_CVX")
cmp


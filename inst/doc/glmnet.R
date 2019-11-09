## ---- eval=FALSE---------------------------------------------------------
#  install.packages("glmnet", repos = "https://cran.us.r-project.org")

## ------------------------------------------------------------------------
library(glmnet)

## ------------------------------------------------------------------------
data(QuickStartExample)

## ------------------------------------------------------------------------
fit = glmnet(x, y)

## ------------------------------------------------------------------------
plot(fit)

## ----height = 4----------------------------------------------------------
print(fit)

## ------------------------------------------------------------------------
coef(fit,s=0.1)

## ------------------------------------------------------------------------
set.seed(29)
nx = matrix(rnorm(10*20),10,20)
predict(fit,newx=nx,s=c(0.1,0.05))

## ------------------------------------------------------------------------
cvfit = cv.glmnet(x, y)

## ------------------------------------------------------------------------
plot(cvfit)

## ------------------------------------------------------------------------
cvfit$lambda.min

## ------------------------------------------------------------------------
coef(cvfit, s = "lambda.min")

## ------------------------------------------------------------------------
predict(cvfit, newx = x[1:5,], s = "lambda.min")

## ------------------------------------------------------------------------
fit = glmnet(x, y, alpha = 0.2, weights = c(rep(1,50),rep(2,50)), nlambda = 20)

## ------------------------------------------------------------------------
print(fit)

## ------------------------------------------------------------------------
plot(fit, xvar = "lambda", label = TRUE)

## ------------------------------------------------------------------------
plot(fit, xvar = "dev", label = TRUE)

## ------------------------------------------------------------------------
fit = glmnet(x, y)
any(fit$lambda == 0.5)
coef.apprx = coef(fit, s = 0.5, exact = FALSE)
coef.exact = coef(fit, s = 0.5, exact = TRUE, x=x, y=y)
cbind2(coef.exact, coef.apprx)

## ------------------------------------------------------------------------
predict(fit, newx = x[1:5,], type = "response", s = 0.05)

## ------------------------------------------------------------------------
cvfit = cv.glmnet(x, y, type.measure = "mse", nfolds = 20)

## ---- eval=FALSE---------------------------------------------------------
#  require(doMC)
#  registerDoMC(cores=2)
#  X = matrix(rnorm(1e4 * 200), 1e4, 200)
#  Y = rnorm(1e4)

## ---- eval=FALSE---------------------------------------------------------
#  system.time(cv.glmnet(X, Y))

## ---- echo=FALSE---------------------------------------------------------
structure(c(2.44, 0.08, 2.518, 0, 0), class = "proc_time", .Names = c("user.self",
"sys.self", "elapsed", "user.child", "sys.child"))

## ---- eval=FALSE---------------------------------------------------------
#  system.time(cv.glmnet(X, Y, parallel = TRUE))

## ---- echo=FALSE---------------------------------------------------------
structure(c(0.508999999999999, 0.057, 1.56699999999999, 1.941,
0.1), class = "proc_time", .Names = c("user.self", "sys.self",
"elapsed", "user.child", "sys.child"))

## ------------------------------------------------------------------------
cvfit$lambda.min
coef(cvfit, s = "lambda.min")
predict(cvfit, newx = x[1:5,], s = "lambda.min")

## ------------------------------------------------------------------------
foldid=sample(1:10,size=length(y),replace=TRUE)
cv1=cv.glmnet(x,y,foldid=foldid,alpha=1)
cv.5=cv.glmnet(x,y,foldid=foldid,alpha=.5)
cv0=cv.glmnet(x,y,foldid=foldid,alpha=0)

## ------------------------------------------------------------------------
par(mfrow=c(2,2))
plot(cv1);plot(cv.5);plot(cv0)
plot(log(cv1$lambda),cv1$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1$name)
points(log(cv.5$lambda),cv.5$cvm,pch=19,col="grey")
points(log(cv0$lambda),cv0$cvm,pch=19,col="blue")
legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("red","grey","blue"))

## ------------------------------------------------------------------------
tfit=glmnet(x,y,lower=-.7,upper=.5)
plot(tfit)

## ------------------------------------------------------------------------
p.fac = rep(1, 20)
p.fac[c(5, 10, 15)] = 0
pfit = glmnet(x, y, penalty.factor = p.fac)
plot(pfit, label = TRUE)

## ------------------------------------------------------------------------
set.seed(101)
x=matrix(rnorm(1000),100,10)
y=rnorm(100)
vn=paste("var",1:10)
fit=glmnet(x,y)
plot(fit)

## ------------------------------------------------------------------------
par(mar=c(4.5,4.5,1,4))
plot(fit)
vnat=coef(fit)
vnat=vnat[-1,ncol(vnat)] # remove the intercept, and get the coefficients at the end of the path
axis(4, at=vnat,line=-.5,label=vn,las=1,tick=FALSE, cex.axis=0.5)

## ------------------------------------------------------------------------
data(MultiGaussianExample)

## ------------------------------------------------------------------------
mfit = glmnet(x, y, family = "mgaussian")

## ------------------------------------------------------------------------
plot(mfit, xvar = "lambda", label = TRUE, type.coef = "2norm")

## ------------------------------------------------------------------------
predict(mfit, newx = x[1:5,], s = c(0.1, 0.01))

## ------------------------------------------------------------------------
cvmfit = cv.glmnet(x, y, family = "mgaussian")

## ------------------------------------------------------------------------
plot(cvmfit)

## ------------------------------------------------------------------------
cvmfit$lambda.min
cvmfit$lambda.1se

## ------------------------------------------------------------------------
data(BinomialExample)

## ------------------------------------------------------------------------
fit = glmnet(x, y, family = "binomial")

## ------------------------------------------------------------------------
plot(fit, xvar = "dev", label = TRUE)

## ------------------------------------------------------------------------
predict(fit, newx = x[1:5,], type = "class", s = c(0.05, 0.01))

## ------------------------------------------------------------------------
cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "class")

## ------------------------------------------------------------------------
plot(cvfit)

## ------------------------------------------------------------------------
cvfit$lambda.min
cvfit$lambda.1se

## ------------------------------------------------------------------------
coef(cvfit, s = "lambda.min")

## ------------------------------------------------------------------------
predict(cvfit, newx = x[1:10,], s = "lambda.min", type = "class")

## ------------------------------------------------------------------------
data(MultinomialExample)

## ------------------------------------------------------------------------
fit = glmnet(x, y, family = "multinomial", type.multinomial = "grouped")

## ------------------------------------------------------------------------
plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")

## ------------------------------------------------------------------------
cvfit=cv.glmnet(x, y, family="multinomial", type.multinomial = "grouped", parallel = TRUE)
plot(cvfit)

## ------------------------------------------------------------------------
predict(cvfit, newx = x[1:10,], s = "lambda.min", type = "class")

## ------------------------------------------------------------------------
data(PoissonExample)

## ------------------------------------------------------------------------
fit = glmnet(x, y, family = "poisson")

## ------------------------------------------------------------------------
plot(fit)

## ------------------------------------------------------------------------
coef(fit, s = 1)
predict(fit, newx = x[1:5,], type = "response", s = c(0.1,1))

## ------------------------------------------------------------------------
cvfit = cv.glmnet(x, y, family = "poisson")

## ------------------------------------------------------------------------
plot(cvfit)

## ------------------------------------------------------------------------
opt.lam = c(cvfit$lambda.min, cvfit$lambda.1se)
coef(cvfit, s = opt.lam)

## ------------------------------------------------------------------------
data(CoxExample)
y[1:5,]

## ------------------------------------------------------------------------
fit = glmnet(x, y, family = "cox")

## ------------------------------------------------------------------------
plot(fit)

## ------------------------------------------------------------------------
coef(fit, s = 0.05)

## ------------------------------------------------------------------------
cvfit = cv.glmnet(x, y, family = "cox")

## ------------------------------------------------------------------------
plot(cvfit)

## ------------------------------------------------------------------------
cvfit$lambda.min
cvfit$lambda.1se

## ------------------------------------------------------------------------
coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
index.min = coef.min[active.min]

## ------------------------------------------------------------------------
index.min
coef.min

## ------------------------------------------------------------------------
data(SparseExample)

## ------------------------------------------------------------------------
class(x)

## ------------------------------------------------------------------------
fit = glmnet(x, y)

## ------------------------------------------------------------------------
cvfit = cv.glmnet(x, y)
plot(cvfit)

## ------------------------------------------------------------------------
i = sample(1:5, size = 25, replace = TRUE)
j = sample(1:20, size = 25, replace = TRUE)
x = rnorm(25)
nx = sparseMatrix(i = i, j = j, x = x, dims = c(5, 20))
predict(cvfit, newx = nx, s = "lambda.min")

## ------------------------------------------------------------------------
data(QuickStartExample)
fit = glmnet(x, y)
print(fit)

## ------------------------------------------------------------------------
glmnet.control(fdev = 0)
fit = glmnet(x, y)
print(fit)

## ------------------------------------------------------------------------
glmnet.control(factory = TRUE)

## ------------------------------------------------------------------------
glmnet.control()

## ---- echo=FALSE---------------------------------------------------------
data(QuickStartExample)

## ----eval=FALSE----------------------------------------------------------
#  fit = glmnet(x, y, intercept = F, standardize = F, lambda = 8/(2*dim(x)[1]), thresh = 1e-20)

## ----eval=FALSE----------------------------------------------------------
#  beta_glmnet = as.matrix(predict(fit, type = "coefficients")[-1,])

## ------------------------------------------------------------------------
fit = glmnet(x, y, intercept = F, standardize = F, thresh = 1e-20)
beta_glmnet = as.matrix(predict(fit, s = 8/(2*dim(x)[1]), type = "coefficients",
                                exact = TRUE, x=x, y=y)[-1,])

## ---- eval=FALSE---------------------------------------------------------
#  library(CVXfromR)
#  setup.dir = "change/this/to/your/cvx/directory"
#  n = dim(x)[1]; p = dim(x)[2]
#  cvxcode = paste("variables beta(p)",
#                  "minimize(square_pos(norm(y - x * beta, 2)) + lambda * norm(beta, 1))",
#                  sep = ";")
#  Lasso = CallCVX(cvxcode, const.var = list(p = p, x = x, y = y, lambda = 8), opt.var.names = "beta", setup.dir = setup.dir, matlab.call = "change/this/to/path/to/matlab")
#  beta_CVX = Lasso$beta

## ------------------------------------------------------------------------
data(CVXResults)

## ----message=FALSE-------------------------------------------------------
require(lars)

## ------------------------------------------------------------------------
fit_lars = lars(x, y, type = "lasso", intercept = F, normalize = F)
beta_lars = predict(fit_lars, s = 8/2, type = "coefficients", mode = "lambda")$coefficients

## ------------------------------------------------------------------------
cmp = round(cbind(beta_glmnet, beta_lars, beta_CVX), digits = 6)
colnames(cmp) = c("beta_glmnet", "beta_lars", "beta_CVX")
cmp


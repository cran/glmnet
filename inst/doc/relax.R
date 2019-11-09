## ------------------------------------------------------------------------
library(glmnet)
data(QuickStartExample)
fit=glmnet(x,y, relax=TRUE)
print(fit)

## ------------------------------------------------------------------------
par(mfrow=c(1,3))
plot(fit)
plot(fit,gamma=0.5)
plot(fit,gamma=0)

## ------------------------------------------------------------------------
cfit=cv.glmnet(x,y,relax=TRUE)
plot(cfit)

## ---- eval=FALSE---------------------------------------------------------
#  predict(cvfit,newx)

## ------------------------------------------------------------------------
print(cfit)

## ----`relaxed`-----------------------------------------------------------
fit=glmnet(x,y)
fitr=relax.glmnet(fit,x=x,y=y)

## ------------------------------------------------------------------------
print(cfit)
print.cv.glmnet(cfit)

## ------------------------------------------------------------------------
fitr=cv.glmnet(x,y,gamma=0,relax=TRUE)
plot(fitr)

## ---- eval=FALSE---------------------------------------------------------
#  fit=glmnet(x,y,trace=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  fit=cv.glmnet(x,y,trace=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  glmnet.control(itrace=1)

## ------------------------------------------------------------------------
 data(CoxExample)

## ------------------------------------------------------------------------
 cvfit=cv.glmnet(x,y,family="cox",type.measure="C")
 plot(cvfit)

## ------------------------------------------------------------------------
data(BinomialExample)
itrain=1:70
fit=glmnet(x[itrain,],y[itrain],family="binomial",nlambda=20)
assess.glmnet(fit,newx=x[-itrain,],newy=y[-itrain])

## ---- eval=FALSE---------------------------------------------------------
#  pred=predict(fit,newx=x[-itrain,])
#  assess.glmnet(pred,newy=y[-itrain],family="binomial")

## ------------------------------------------------------------------------
glmnet.measures()

## ------------------------------------------------------------------------
cfit=cv.glmnet(x[itrain,],y[itrain],family="binomial", nlambda = 30)
assess.glmnet(cfit,newx=x[-itrain,],newy=y[-itrain])

## ------------------------------------------------------------------------
assess.glmnet(cfit,newx=x[-itrain,],newy=y[-itrain], s="lambda.min")

## ------------------------------------------------------------------------
cfit=cv.glmnet(x,y,family="binomial",keep=TRUE, nlambda = 30)
assess.glmnet(cfit$fit.preval,newy=y,family="binomial")

## ------------------------------------------------------------------------
cfit=cv.glmnet(x,y,family="binomial", type.measure="auc", keep=TRUE)
rocs=roc.glmnet(cfit$fit.preval,newy=y)
which=match(cfit$lambda.min,cfit$lambda)
plot(rocs[[which]],type="l")
nopr=sapply(rocs,lines,col="grey")
lines(rocs[[which]],lwd=2,col="red")

## ------------------------------------------------------------------------
data(MultinomialExample)
set.seed(101)
itrain=sample(1:500,400,replace=FALSE)
cfit=cv.glmnet(x[itrain,],y[itrain],family="multinomial")
cnf=confusion.glmnet(cfit,newx=x[-itrain,],newy=y[-itrain])
print(cnf)

## ------------------------------------------------------------------------
cfit=cv.glmnet(x,y,family="multinomial",type="class",keep=TRUE)
cnf=confusion.glmnet(cfit$fit.preval,newy=y,family="multinomial")
which=match(cfit$lambda.min,cfit$lambda)
print(cnf[[which]])

## ------------------------------------------------------------------------
data(BinomialExample)
fit=bigGlm(x,y,family="binomial",lower.limits=-1)
print(fit)

## ------------------------------------------------------------------------
set.seed(101)
X = matrix(rnorm(20),10,2)
X3=sample(letters[1:3],10,replace=TRUE)
X4=sample(LETTERS[1:3],10,replace=TRUE)
df=data.frame(X,X3,X4)
makeX(df)

## ------------------------------------------------------------------------
makeX(df,sparse=TRUE)

## ------------------------------------------------------------------------
Xn=X
Xn[3,1]=NA;Xn[5,2]=NA
X3n=X3;
X3n[6]=NA
X4n=X4
X4n[9]=NA
dfn=data.frame(Xn,X3n,X4n)
makeX(dfn)

## ------------------------------------------------------------------------
makeX(dfn,na.impute=TRUE,sparse=TRUE)

## ------------------------------------------------------------------------
X = matrix(rnorm(10),5,2)
X3=sample(letters[1:3],5,replace=TRUE)
X4=sample(LETTERS[1:3],5,replace=TRUE)
Xn=X
Xn[3,1]=NA;Xn[5,2]=NA
X3n=X3;
X3n[1]=NA
X4n=X4
X4n[2]=NA
dftn=data.frame(Xn,X3n,X4n)
makeX(dfn,dftn,na.impute=TRUE, sparse=TRUE)


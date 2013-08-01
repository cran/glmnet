### R code from vignette source 'Coxnet.rnw'

###################################################
### code chunk number 1: Coxnet.rnw:30-34
###################################################
library("glmnet")
library("survival")
load("VignetteExample.rdata")
attach(patient.data)


###################################################
### code chunk number 2: Coxnet.rnw:44-46
###################################################
cv.fit <- cv.glmnet(x, Surv(time, status), family="cox", maxit = 1000)
fit <- glmnet(x, Surv(time,status), family =  "cox", maxit = 1000)


###################################################
### code chunk number 3: Coxnet.rnw:56-57
###################################################
plot(cv.fit)


###################################################
### code chunk number 4: Coxnet.rnw:60-61
###################################################
cv.fit$lambda.min


###################################################
### code chunk number 5: Coxnet.rnw:74-77
###################################################
Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients  <- Coefficients[Active.Index]


###################################################
### code chunk number 6: Coxnet.rnw:82-84
###################################################
Active.Index
Active.Coefficients



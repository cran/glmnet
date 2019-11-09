## ------------------------------------------------------------------------
library("glmnet")
library("survival")
patient.data  <- readRDS("assets/coxnet.RDS")

## ---- warning = TRUE-----------------------------------------------------
cv.fit <- cv.glmnet(patient.data$x, Surv(patient.data$time, patient.data$status), family="cox", maxit = 1000)
fit <- glmnet(patient.data$x, Surv(patient.data$time,patient.data$status), family =  "cox", maxit = 1000)

## ------------------------------------------------------------------------
plot(cv.fit)
cv.fit$lambda.min

## ------------------------------------------------------------------------
Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients  <- Coefficients[Active.Index]

## ------------------------------------------------------------------------
Active.Index
Active.Coefficients


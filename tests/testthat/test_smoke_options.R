## CRAN smoke tests: exercise each non-trivial glmnet() option once.
## Uses the smaller (n=500, p=50) setup from test_smoke_families.R.
## Fixture-free; reproducible via set.seed.

library(glmnet)

SEED <- 20260421L
N <- 500L
P <-  50L

set.seed(SEED)
X <- matrix(rnorm(N * P), N, P)
set.seed(SEED + 100L)
y_gauss <- rnorm(N)
set.seed(SEED + 101L)
y_bin   <- rbinom(N, 1, 0.5)

test_that("alpha grid (gaussian): 0, 0.5, 1", {
  for (a in c(0, 0.5, 1)) {
    fit <- glmnet(X, y_gauss, alpha = a)
    expect_s3_class(fit, "glmnet")
    expect_true(length(fit$lambda) >= 1L)
  }
})

test_that("intercept=FALSE fits with a0 == 0 (gaussian)", {
  fit <- glmnet(X, y_gauss, intercept = FALSE)
  expect_true(all(fit$a0 == 0))
})

test_that("observation weights (gaussian)", {
  set.seed(SEED + 200L)
  w <- runif(N, 0.5, 1.5)
  fit <- glmnet(X, y_gauss, weights = w)
  expect_s3_class(fit, "glmnet")
  expect_equal(ncol(fit$beta), length(fit$lambda))
})

test_that("offset (poisson)", {
  set.seed(SEED + 300L)
  y_pois <- rpois(N, 2)
  off <- runif(N, -0.2, 0.2)
  fit <- glmnet(X, y_pois, family = "poisson", offset = off)
  expect_s3_class(fit, "glmnet")
  pr <- predict(fit, X, s = min(fit$lambda), newoffset = off,
                type = "response")
  expect_true(all(pr > 0))
})

test_that("exclude excludes variables from the model", {
  ex <- c(1L, 5L, 20L)
  fit <- glmnet(X, y_gauss, exclude = ex)
  expect_true(all(fit$beta[ex, ] == 0))
})

test_that("penalty.factor with zero for unpenalized variables", {
  pf <- rep(1, P); pf[c(3L, 7L)] <- 0
  fit <- glmnet(X, y_gauss, penalty.factor = pf)
  expect_s3_class(fit, "glmnet")
  # unpenalized variables should be active at the largest lambda
  expect_true(any(fit$beta[3L, ] != 0))
  expect_true(any(fit$beta[7L, ] != 0))
})

test_that("lower.limits / upper.limits clip coefficients", {
  lb <- rep(-0.1, P); ub <- rep(0.1, P)
  fit <- glmnet(X, y_gauss, lower.limits = lb, upper.limits = ub)
  expect_true(all(as.matrix(fit$beta) >= -0.1 - 1e-8))
  expect_true(all(as.matrix(fit$beta) <=  0.1 + 1e-8))
})

test_that("family = gaussian() goes through glmnet.path / glmnet.fit", {
  fit <- glmnet(X, y_gauss, family = gaussian())
  expect_s3_class(fit, c("glmnetfit", "glmnet"))
  expect_equal(ncol(fit$beta), length(fit$lambda))
})

test_that("family = binomial() goes through glmnet.path / glmnet.fit", {
  fit <- glmnet(X, y_bin, family = binomial())
  expect_s3_class(fit, c("glmnetfit", "glmnet"))
  prob <- predict(fit, X, s = min(fit$lambda), type = "response")
  expect_true(all(prob >= 0 & prob <= 1))
})

test_that("control = list() per-call override of thresh", {
  fit_tight <- glmnet(X, y_gauss, control = list(thresh = 1e-10))
  fit_loose <- glmnet(X, y_gauss, control = list(thresh = 1e-3))
  expect_s3_class(fit_tight, "glmnet")
  expect_s3_class(fit_loose, "glmnet")
})

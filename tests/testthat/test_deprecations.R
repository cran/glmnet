## CRAN smoke tests: deprecation + transition warnings.

library(glmnet)

SEED <- 20260421L
N <- 100L; P <- 20L
set.seed(SEED)
X <- matrix(rnorm(N * P), N, P)

test_that("cox.ties unset -> transition warning (v5.0 -> v5.1)", {
  set.seed(SEED + 1L)
  time   <- runif(N, 1, 10)
  status <- rbinom(N, 1, 0.7)
  y <- survival::Surv(time, status)
  expect_warning(glmnet(X, y, family = "cox"),
                 regexp = "cox\\.ties")
  ## Explicitly setting cox.ties suppresses the warning.
  expect_silent(glmnet(X, y, family = "cox", cox.ties = "breslow"))
  expect_silent(glmnet(X, y, family = "cox", cox.ties = "efron"))
})

test_that("top-level thresh/maxit/dfmax/pmax/trace.it emit .Deprecated", {
  set.seed(SEED + 2L)
  y <- rnorm(N)
  expect_warning(glmnet(X, y, thresh   = 1e-6),  regexp = "deprecated")
  expect_warning(glmnet(X, y, maxit    = 1000L), regexp = "deprecated")
  expect_warning(glmnet(X, y, dfmax    = 10L),   regexp = "deprecated")
  expect_warning(glmnet(X, y, pmax     = 15L),   regexp = "deprecated")
  expect_warning(glmnet(X, y, trace.it = 0L),    regexp = "deprecated")
})

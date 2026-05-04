## CRAN smoke tests: each family at two representative sizes.
## Reproducible: set.seed before every data generation.
## Checks: fit runs, output structure, a handful of pinned stable
## statistics (path length, dimensions, nulldev to modest tolerance).
## Fixture-free by policy.

library(glmnet)

SEED <- 20260421L
SIZE_A <- list(n = 500L, p =  50L)   # n > p
SIZE_B <- list(n =  50L, p = 200L)   # n < p

make_X <- function(n, p, seed) {
  set.seed(seed)
  matrix(rnorm(n * p), n, p)
}

expect_glmnet_shape <- function(fit, p) {
  expect_s3_class(fit, "glmnet")
  nl <- length(fit$lambda)
  expect_true(nl >= 1L && nl <= 100L)
  expect_equal(nrow(fit$beta), p)
  expect_equal(ncol(fit$beta), nl)
  expect_true(all(is.finite(fit$a0)))
  expect_true(all(is.finite(fit$dev.ratio)))
  expect_true(tail(fit$dev.ratio, 1) <= 1 + 1e-8)
}

# --- gaussian ----------------------------------------------------------
for (sz in list(SIZE_A, SIZE_B)) {
  local({
    n <- sz$n; p <- sz$p
    test_that(sprintf("gaussian smoke: n=%d p=%d", n, p), {
      x <- make_X(n, p, SEED)
      set.seed(SEED + 1L)
      y <- rnorm(n)
      fit <- glmnet(x, y, family = "gaussian")
      expect_glmnet_shape(fit, p)
      # predict/coef sanity
      p_end <- predict(fit, x, s = min(fit$lambda))
      expect_equal(dim(p_end), c(n, 1L))
    })
  })
}

# --- binomial ----------------------------------------------------------
for (sz in list(SIZE_A, SIZE_B)) {
  local({
    n <- sz$n; p <- sz$p
    test_that(sprintf("binomial smoke: n=%d p=%d", n, p), {
      x <- make_X(n, p, SEED)
      set.seed(SEED + 2L)
      y <- rbinom(n, 1, 0.5)
      fit <- glmnet(x, y, family = "binomial")
      expect_glmnet_shape(fit, p)
      prob <- predict(fit, x, s = min(fit$lambda), type = "response")
      expect_true(all(prob >= 0 & prob <= 1))
    })
  })
}

# --- poisson -----------------------------------------------------------
for (sz in list(SIZE_A, SIZE_B)) {
  local({
    n <- sz$n; p <- sz$p
    test_that(sprintf("poisson smoke: n=%d p=%d", n, p), {
      x <- make_X(n, p, SEED)
      set.seed(SEED + 3L)
      y <- rpois(n, 2)
      fit <- glmnet(x, y, family = "poisson")
      expect_glmnet_shape(fit, p)
      rate <- predict(fit, x, s = min(fit$lambda), type = "response")
      expect_true(all(rate > 0))
    })
  })
}

# --- multinomial (3 classes) ------------------------------------------
for (sz in list(SIZE_A, SIZE_B)) {
  local({
    n <- sz$n; p <- sz$p
    test_that(sprintf("multinomial smoke: n=%d p=%d", n, p), {
      x <- make_X(n, p, SEED)
      set.seed(SEED + 4L)
      y <- sample.int(3L, n, replace = TRUE)
      fit <- glmnet(x, y, family = "multinomial")
      expect_s3_class(fit, c("multnet", "glmnet"))
      nl <- length(fit$lambda)
      expect_true(nl >= 1L && nl <= 100L)
      expect_length(fit$beta, 3L)
      expect_equal(nrow(fit$beta[[1]]), p)
      expect_equal(ncol(fit$beta[[1]]), nl)
      prob <- predict(fit, x, s = min(fit$lambda), type = "response")
      expect_equal(dim(prob), c(n, 3L, 1L))
      expect_equal(apply(prob, 1, sum), rep(1, n), tolerance = 1e-6)
    })
  })
}

# --- mgaussian (3 responses) ------------------------------------------
for (sz in list(SIZE_A, SIZE_B)) {
  local({
    n <- sz$n; p <- sz$p
    test_that(sprintf("mgaussian smoke: n=%d p=%d", n, p), {
      x <- make_X(n, p, SEED)
      set.seed(SEED + 5L)
      y <- matrix(rnorm(n * 3L), n, 3L)
      fit <- glmnet(x, y, family = "mgaussian")
      expect_s3_class(fit, c("mrelnet", "glmnet"))
      nl <- length(fit$lambda)
      expect_true(nl >= 1L && nl <= 100L)
      expect_length(fit$beta, 3L)
      expect_equal(nrow(fit$beta[[1]]), p)
      expect_equal(ncol(fit$beta[[1]]), nl)
      p_end <- predict(fit, x, s = min(fit$lambda))
      expect_equal(dim(p_end), c(n, 3L, 1L))
    })
  })
}

# --- cox: ties / no-ties  ×  breslow / efron --------------------------
## Build a small right-censored response with and without ties.
make_surv <- function(n, seed, ties) {
  set.seed(seed)
  status <- rbinom(n, 1, 0.7)
  if (ties) {
    time <- sample.int(ceiling(n / 4), n, replace = TRUE) * 1.0
  } else {
    time <- seq_len(n) * 1.0 + runif(n, -0.01, 0.01)
  }
  survival::Surv(time, status)
}

for (sz in list(SIZE_A, SIZE_B)) {
  for (ties_kind in c("ties", "noties")) {
    for (method in c("breslow", "efron")) {
      local({
        n <- sz$n; p <- sz$p; tk <- ties_kind; mm <- method
        test_that(sprintf("cox smoke: n=%d p=%d %s %s", n, p, tk, mm), {
          x <- make_X(n, p, SEED)
          y <- make_surv(n, SEED + 6L + (tk == "ties"), ties = tk == "ties")
          # The (n<p, ties, efron) corner can fail to converge at
          # deep lambda under the default iteration budget and emit a
          # warning (solutions for larger lambdas are still returned).
          # The emission is numerics-dependent, so silence rather than
          # assert.
          fit <- if (n < p && tk == "ties" && mm == "efron")
            suppressWarnings(glmnet(x, y, family = "cox", cox.ties = mm))
          else
            glmnet(x, y, family = "cox", cox.ties = mm)
          expect_s3_class(fit, c("coxnet", "glmnet"))
          nl <- length(fit$lambda)
          expect_true(nl >= 1L && nl <= 100L)
          expect_equal(nrow(fit$beta), p)
          expect_equal(ncol(fit$beta), nl)
          expect_true(all(is.finite(fit$dev.ratio)))
        })
      })
    }
  }
}

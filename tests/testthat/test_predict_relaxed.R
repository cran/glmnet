## predict() / coef() on relaxed fits with a vector of gamma values must
## return a named list equal to the corresponding single-gamma calls.
## (This path was broken from its introduction until 5.2.)

library(glmnet)

set.seed(20261020L)
N <- 200L
P <-  20L
X <- matrix(rnorm(N * P), N, P)
eta <- drop(X[, 1:4] %*% c(1, -1, 0.5, 0.5))
y_gauss <- eta + rnorm(N)
y_bin   <- rbinom(N, 1, plogis(eta))
y_mult  <- factor(apply(cbind(eta, -eta, 0) + matrix(rexp(3 * N), N), 1, which.max))
off     <- rnorm(N, sd = 0.2)
newx    <- X[1:5, ]
gammas  <- c(0, 0.25, 1)

expect_matches_single <- function(multi, single_fun) {
  expect_type(multi, "list")
  expect_named(multi, format(round(gammas, 2)))
  for (i in seq_along(gammas))
    expect_equal(multi[[i]], single_fun(gammas[i]))
}

test_that("vector gamma: gaussian link predictions and coefficients", {
  fit <- glmnet(X, y_gauss, relax = TRUE)
  expect_matches_single(
    predict(fit, newx = newx, s = c(0.05, 0.1), gamma = gammas),
    function(g) predict(fit, newx = newx, s = c(0.05, 0.1), gamma = g))
  expect_matches_single(
    coef(fit, s = 0.1, gamma = gammas),
    function(g) coef(fit, s = 0.1, gamma = g))
  expect_matches_single(
    predict(fit, s = 0.1, gamma = gammas, type = "nonzero"),
    function(g) predict(fit, s = 0.1, gamma = g, type = "nonzero"))
})

test_that("vector gamma: binomial type is passed through", {
  fit <- glmnet(X, y_bin, family = "binomial", relax = TRUE)
  for (ty in c("link", "response", "class"))
    expect_matches_single(
      predict(fit, newx = newx, s = 0.02, gamma = gammas, type = ty),
      function(g) predict(fit, newx = newx, s = 0.02, gamma = g, type = ty))
  ## response must differ from link, i.e. type was not dropped
  r <- predict(fit, newx = newx, s = 0.02, gamma = gammas, type = "response")
  expect_true(all(r[[1]] > 0 & r[[1]] < 1))
})

test_that("vector gamma: multinomial (list of coefficient matrices)", {
  fit <- glmnet(X, y_mult, family = "multinomial", relax = TRUE)
  expect_matches_single(
    predict(fit, newx = newx, s = 0.02, gamma = gammas, type = "class"),
    function(g) predict(fit, newx = newx, s = 0.02, gamma = g, type = "class"))
  expect_matches_single(
    coef(fit, s = 0.02, gamma = gammas),
    function(g) coef(fit, s = 0.02, gamma = g))
})

test_that("vector gamma: newoffset is passed through for offset fits", {
  fit <- glmnet(X, y_bin, family = "binomial", offset = off, relax = TRUE)
  expect_matches_single(
    predict(fit, newx = newx, s = 0.02, gamma = gammas, newoffset = off[1:5]),
    function(g) predict(fit, newx = newx, s = 0.02, gamma = g, newoffset = off[1:5]))
})

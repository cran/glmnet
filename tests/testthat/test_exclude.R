# TESTS FOR EXCLUDE ARGUMENT

# passes when everything matches
compare_path_fits <- function(oldfit, newfit) {
    expect_equal(coef(oldfit), coef(newfit), tolerance = 1e-07,
                 label = "coef doesn't match")
    expect_equal(oldfit$df, newfit$df,
                 label = "df doesn't match")
    expect_equal(oldfit$dim, newfit$dim,
                 label = "dim doesn't match")
    expect_equal(oldfit$lambda, newfit$lambda,
                 label = "lambda doesn't match")
    expect_equal(oldfit$dev.ratio, newfit$dev.ratio,
                 label = "dev.ratio doesn't match")
    expect_equal(oldfit$nulldev, newfit$nulldev,
                 label = "nulldev doesn't match")
    expect_equal(oldfit$offset, newfit$offset,
                 label = "offset flag doesn't match")
    expect_equal(oldfit$nobs, newfit$nobs,
                 label = "nobs doesn't match")
}


# parameters for data (and weights for glmnet & data)
nobs <- 100; nvars <- 15
beta <- matrix(c(2, -2, rep(0, nvars - 2)), ncol = 1)

# set up fake data
set.seed(4)
x <- matrix(rnorm(nobs * nvars), nrow = nobs)
y <- (x %*% beta + rnorm(nobs)) / 3
weights <- rep(1:2, length.out = nobs)

test_that("No variables excluded", {
    exclude_fit <- glmnet(x, y, exclude = c())
    glmnet_fit <- glmnet(x, y)
    compare_path_fits(exclude_fit, glmnet_fit)
})

test_that("Exclude defined by fixed vector", {
  exclude_fit <- glmnet(x, y, exclude = 1:5)
  
  x2 <- x[, 6:nvars]
  glmnet_fit <- glmnet(x2, y)
  glmnet_fit$beta <- rbind(matrix(0, nrow = 5, ncol = ncol(glmnet_fit$beta)),
                           glmnet_fit$beta)
  rownames(glmnet_fit$beta) <- paste0("V", 1:nvars)
  glmnet_fit$dim <- c(nvars, length(glmnet_fit$lambda))
  compare_path_fits(exclude_fit, glmnet_fit)
})

test_that("Exclude defined by function of x, no variables filtered", {
  filter_fn <- function(x, ...) which(colMeans(x == 0) > 0.8)
  exclude_fit <- glmnet(x, y, exclude = filter_fn)
  glmnet_fit <- glmnet(x, y)
  compare_path_fits(exclude_fit, glmnet_fit)
})

test_that("Exclude defined by function of x, some variables filtered", {
  # for our test data, this excludes columns 3, 5, 6, 9, 11
  filter_fn <- function(x, ...) which(apply(x, 2, max) > 2.5)
  exclude_fit <- glmnet(x, y, exclude = filter_fn)
  
  exclude <- c(3, 5, 6, 9, 11)
  exclude_fit2 <- glmnet(x, y, exclude = exclude)
  compare_path_fits(exclude_fit, exclude_fit2)
  
  include <- setdiff(1:nvars, exclude)
  x2 <- x[, include]
  glmnet_fit <- glmnet(x2, y)
  temp_beta <- Matrix::Matrix(0, nrow = nvars, ncol = ncol(glmnet_fit$beta),
                              sparse = TRUE)
  temp_beta[include, ] <- glmnet_fit$beta
  rownames(temp_beta) <- paste0("V", 1:nvars)
  glmnet_fit$beta <- temp_beta
  glmnet_fit$dim <- c(nvars, length(glmnet_fit$lambda))
  compare_path_fits(exclude_fit, glmnet_fit)
})

test_that("Exclude defined by function of x and y, some variables filtered", {
  # for our test data, this excludes all columns except 1 and 2
  filter_fn <- function(x, y, ...) which(apply(x, 2, 
                                               function(x) abs(cor(x, y))) <= 0.5)
  exclude_fit <- glmnet(x, y, exclude = filter_fn)
  
  exclude <- 3:nvars
  exclude_fit2 <- glmnet(x, y, exclude = exclude)
  compare_path_fits(exclude_fit, exclude_fit2)
})

test_that("Exclude defined by function of x and weights, some variables filtered", {
  # for our test data, this excludes all columns except 2, 5 and 13
  filter_fn <- function(x, weights, ...) which(
    apply(x, 2, function(x) weighted.mean((x - mean(x))^2, weights)) <= 1)
  exclude_fit <- glmnet(x, y, weights = weights, exclude = filter_fn)
  
  exclude <- setdiff(1:nvars, c(2, 5, 13))
  exclude_fit2 <- glmnet(x, y, weights = weights, exclude = exclude)
  compare_path_fits(exclude_fit, exclude_fit2)
})

test_that("cv.glmnet with exclude defined by function of x", {
  # for this, we only compare the fit.preval values
  filter_fn <- function(x, y, ...) which(apply(x, 2, 
                                               function(x) abs(cor(x, y))) <= 0.1)
  foldid <- rep(1:4, each = 25)
  
  exclude_fit <- cv.glmnet(x, y, exclude = filter_fn, foldid = foldid, 
                           keep = TRUE)
  
  # do CV by hand
  overall_exclude <- filter_fn(x, y)
  glmnet_fit <- glmnet(x, y, exclude = overall_exclude)
  fit.preval_by_hand <- matrix(NA, nrow = nobs, ncol = length(glmnet_fit$lambda))
  colnames(fit.preval_by_hand) <- paste0("s", 1:ncol(fit.preval_by_hand) - 1)
  for (k in 1:4) {
    in_idx <- which(foldid != k)
    out_idx <- which(foldid == k)
    exclude <- filter_fn(x[in_idx, ], y[in_idx])
    fit <- glmnet(x[in_idx, ], y[in_idx], exclude = exclude)
    fit.preval_by_hand[out_idx, ] <- predict(fit, newx = x[out_idx, ], 
                                             s = glmnet_fit$lambda)
  }
  expect_equal(exclude_fit$fit.preval, fit.preval_by_hand)
})

#####
# in the future, these calls may not throw errors
#####
test_that("All variables excluded", {

  expect_error(glmnet(x, y, exclude = 1:nvars))
})

# Test for glmnet.path with new families
#
# All tests are performed for intercept = TRUE and FALSE and standardize = TRUE
# and FALSE (except subgradient equations, which are only performed for
# intercept = TRUE and FALSE).
#
# binomial with probit link
# binomial with probit link & observation weights
# binomial with probit link, subgradient equations with obs weights
# quasipoisson & observation weights
# quasipoisson subgradient equations with obs weights

compare_fits <- function(oldfit, newfit, intr) {
    if (intr) {
        expect_equal(as.numeric(coef(oldfit)), as.numeric(coef(newfit)),
                     label = "coefficients don't match")
    } else {
        expect_equal(as.numeric(coef(oldfit)), as.numeric(coef(newfit))[-1],
                     label = "coefficients don't match")
    }
    expect_equal(oldfit$null.deviance, newfit$nulldev,
                 label = "nulldev doesn't match")
    expect_equal(1 - oldfit$deviance / oldfit$null.deviance, newfit$dev.ratio,
                 label = "dev.ratio doesn't match")
}

check_subgrad <- function(newfit, x, y, weights) {
    nobs <- newfit$nobs; nvars <- nrow(newfit$beta)
    nlam <- length(newfit$lambda)
    weights <- weights / sum(weights)

    mu <- predict(newfit, x, type = "response")
    eta <- predict(newfit, x, type = "link")
    r <- matrix(y, nrow = nobs, ncol = nlam) - mu
    v <- newfit$family$variance(mu)
    m.e <- newfit$family$mu.eta(eta)
    subgrad <- abs(t(r * m.e / v) %*% (x * weights))  # this is nlam by nvars matrix

    for (i in 1:nlam) {
        for (j in 1:nvars) {
            if (abs(newfit$beta[j, i]) > 0) {
                expect_equal(as.numeric(subgrad[i, j]),
                             as.numeric(newfit$lambda[i]),
                             label = paste("lambda", i, "var", j))
            } else {
                # + 1e-15 to allow for a little slack
                expect_lte(subgrad[i, j], newfit$lambda[i] + 1e-15,
                           label = paste("lambda", i, "var", j))
            }
        }
    }
}


# parameters for data
nobs <- 150; nvars <- 10
beta <- matrix(c(1:5, -3, rep(0, nvars - 6)), ncol = 1)

# set up fake data
set.seed(400)
x <- matrix(rnorm(nobs * nvars), nrow = nobs)
y <- (x %*% beta + 5 * rnorm(nobs)) / 5
biny <- ifelse(y > 0, 1, 0)
poiy <- rpois(length(y), exp(y))

# probit link
eps <- 1e-14
thresh <- 1e-16
family <- binomial(link = "probit")
glmnet.control(epsnr = eps)
for (intr in c(TRUE, FALSE)) {
    for (isd in c(TRUE, FALSE)) {
        test_that(paste("binomial probit link, intercept", intr, "std", isd), {
            if (intr) {
                glmfit <- glm(biny ~ x, family = family, control = list(epsilon = eps))
            } else {
                glmfit <- glm(biny ~ x - 1, family = family, control = list(epsilon = eps))
            }
            newfit <- glmnet(x, biny, lambda = 0, intercept = intr, 
                             standardize = isd, family = family, control = .control_with(thresh = thresh))
            compare_fits(glmfit, newfit, intr)
        })
    }
}

# probit link with weights
weights <- rep(1:3, length.out = nobs)
eps <- 1e-14
thresh <- 1e-16
family <- binomial(link = "probit")
glmnet.control(epsnr = eps)
for (intr in c(TRUE, FALSE)) {
    for (isd in c(TRUE, FALSE)) {
        test_that(paste("binomial probit link & obs weights, intercept", intr,
                        "std", isd), {
            if (intr) {
                glmfit <- glm(biny ~ x, family = family, control = list(epsilon = eps),
                              weights = weights)
            } else {
                glmfit <- glm(biny ~ x - 1, family = family, control = list(epsilon = eps),
                              weights = weights)
            }
            newfit <- glmnet(x, biny, lambda = 0, intercept = intr,
                             standardize = FALSE, weights = weights,
                             family = family, control = .control_with(thresh = thresh))
            compare_fits(glmfit, newfit, intr)
        })
    }
}

# probit link: subgradient equations with weights
family <- binomial(link = "probit")
eps <- 1e-14
thresh <- 1e-16
weights <- rep(1:3, length.out = nobs)
glmnet.control(epsnr = eps)
for (intr in c(TRUE, FALSE)) {
    test_that(paste("binomial probit subgrad, intercept", intr), {
        newfit <- glmnet(x, biny, weights = weights,
                         intercept = intr, standardize = FALSE,
                         family = family, control = .control_with(thresh = thresh))
        check_subgrad(newfit, x, biny, weights)
    })
}

# quasipoisson with observation weights
eps <- 1e-14
thresh <- 1e-16
family <- quasipoisson()
weights <- rep(2:5, length.out = nobs)
glmnet.control(epsnr = eps)
for (intr in c(TRUE, FALSE)) {
    for (isd in c(TRUE, FALSE)) {
        test_that(paste("quasipoisson, intercept", intr, "std", isd), {
            if (intr) {
                glmfit <- glm(poiy ~ x, family = family, control = list(epsilon = eps),
                              weights = weights)
            } else {
                glmfit <- glm(poiy ~ x - 1, family = family, control = list(epsilon = eps),
                              weights = weights)
            }
            newfit <- glmnet(x, poiy, lambda = 0, intercept = intr,
                             standardize = isd, weights = weights,
                             family = family, control = .control_with(thresh = thresh))
            compare_fits(glmfit, newfit, intr)
        })
    }
}

# quasipoisson: subgradient equations with weights
family <- quasipoisson()
eps <- 1e-14
thresh <- 1e-16
weights <- rep(1:1, length.out = nobs)
glmnet.control(epsnr = eps)
for (intr in c(TRUE, FALSE)) {
    test_that(paste("quasipoisson subgrad, intercept", intr), {
        newfit <- glmnet(x, poiy, weights = weights,
                         intercept = intr, standardize = FALSE,
                         family = family, control = .control_with(thresh = thresh))
        check_subgrad(newfit, x, poiy, weights)
    })
}

glmnet.control(factory = TRUE)
# TESTS FOR STRATIFYSURV CLASS


y <-Surv(1:10, event = rep(0:1, 5))
strata <- rep(1:2, 5)
y2 <- stratifySurv(y, strata)

expect_equal(class(y2), c("stratifySurv", "Surv"))
expect_equal(attr(y2, "strata"), strata)

y3 <- y2[1]
expect_equal(class(y3), c("stratifySurv", "Surv"))
expect_equal(attr(y3, "strata"), strata[1])

y3 <- y2[2:4]
expect_equal(class(y3), c("stratifySurv", "Surv"))
expect_equal(attr(y3, "strata"), strata[2:4])

y3 <- y2[1, ]
expect_equal(class(y3), c("stratifySurv", "Surv"))
expect_equal(attr(y3, "strata"), strata[1])

y3 <- y2[2:4, ]
expect_equal(class(y3), c("stratifySurv", "Surv"))
expect_equal(attr(y3, "strata"), strata[2:4])

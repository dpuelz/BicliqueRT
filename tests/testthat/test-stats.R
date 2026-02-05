test_that("one_sided_test returns valid p-value", {
  tobs <- 2.0
  tvals <- rnorm(100)
  pval <- one_sided_test(tobs, tvals)
  expect_true(pval >= 0 && pval <= 1)
})

test_that("two_sided_test returns valid p-value", {
  tobs <- 0.5
  tvals <- rnorm(100)
  pval <- two_sided_test(tobs, tvals)
  expect_true(pval >= 0 && pval <= 1)
})

test_that("out_pval returns p-value or decision", {
  rtest_out <- list(tobs = 1.5, tvals = rnorm(100))

  pval <- out_pval(rtest_out, ret_pval = TRUE, alpha = 0.05)
  expect_true(is.numeric(pval))
  expect_true(pval >= 0 && pval <= 2)
})

test_that("gen_tstat returns a function for each type", {
  G <- matrix(0, 5, 5)
  G[1, 2] <- G[2, 1] <- 1
  G[2, 3] <- G[3, 2] <- 1

  for (type in c("elc", "score", "htn")) {
    f <- gen_tstat(G, type)
    expect_true(is.function(f))

    # test that it can be called
    y <- rnorm(5)
    z <- c(1, 0, 1, 0, 1)
    is_focal <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
    val <- f(y, z, is_focal)
    expect_true(is.numeric(val))
  }
})

test_that("ate computes difference in means", {
  Z <- c(1, 1, -1, -1)
  Y <- c(3, 5, 1, 2)
  result <- ate(Z, Y)
  expect_equal(result, mean(c(3, 5)) - mean(c(1, 2)))
})

test_that("out_Yobs computes observed outcomes correctly", {
  Za <- c(1, 1, 0)
  Zb <- c(0, 0, 1)
  Ya <- c(0, 0, 0)
  Yb <- c(1, 1, 1)
  result <- out_Yobs(Za, Zb, Ya, Yb)
  expect_equal(result, c(0, 0, 1))
})

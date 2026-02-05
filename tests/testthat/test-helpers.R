test_that("out_example_network returns correct structure", {
  set.seed(1)
  net <- out_example_network(100)
  expect_true(is.list(net))
  expect_true("D" %in% names(net))
  expect_true("thepoints" %in% names(net))
  expect_equal(nrow(net$D), 100)
  expect_equal(ncol(net$D), 100)
  expect_equal(nrow(net$thepoints), 100)
})

test_that("out_Z generates correct dimensions", {
  pi_vec <- rep(0.5, 10)
  Z <- out_Z(pi_vec, 50)
  expect_equal(nrow(Z), 10)
  expect_equal(ncol(Z), 50)
})

test_that("out_house_structure returns correct structure", {
  hs <- out_house_structure(N = 90, K = 30, equal = TRUE)
  expect_equal(length(hs), 30)
  expect_equal(sum(hs), 90)

  set.seed(1)
  hs2 <- out_house_structure(N = 90, K = 30, equal = FALSE)
  expect_equal(length(hs2), 30)
})

test_that("out_treat_household selects households", {
  hs <- rep(3, 30)
  treat <- out_treat_household(hs, K1 = 15)
  expect_equal(length(treat$which_house), 15)
  expect_equal(length(treat$which_kid), 15)
})

test_that("out_Z_household returns correct structure", {
  N <- 90
  K <- 30
  hs <- rep(3, K)
  treat <- out_treat_household(hs, K1 = 15)
  Z <- out_Z_household(N, K, treat, hs)
  expect_equal(nrow(Z), N)
  expect_equal(ncol(Z), 3)
})

test_that("sparsify returns a sparse matrix", {
  m <- matrix(c(1, 0, 0, 1), nrow = 2)
  sm <- sparsify(m)
  expect_true(inherits(sm, "Matrix") || inherits(sm, "sparseMatrix"))
})

test_that("sample_mod returns integer in range", {
  set.seed(1)
  x <- sample_mod(10)
  expect_true(x >= 1 && x <= 10)
})

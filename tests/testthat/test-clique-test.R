test_that("clique_test runs and returns expected structure", {
  set.seed(123)
  num_units <- 20
  Z <- rbinom(num_units, 1, 0.5)
  Y <- rnorm(num_units)

  hypothesis <- list(
    design_fn = function() rbinom(num_units, 1, 0.5),
    exposure_i = function(z, i) z[i],
    null_equiv = function(e1, e2) identical(e1, e2)
  )
  controls <- list(method = "greedy", mina = 5, num_randomizations = 50)

  # Build adjacency matrix for test stat
  G <- matrix(0, num_units, num_units)
  for (i in 1:(num_units - 1)) G[i, i + 1] <- G[i + 1, i] <- 1

  decomp <- suppressMessages(biclique.decompose(Z, hypothesis, controls))
  tstat <- gen_tstat(G, "htn")
  result <- clique_test(Y, Z, tstat, decomp, alpha = 0.05, one_sided = TRUE)

  expect_true(is.list(result))
  expect_true("p.value" %in% names(result))
  expect_true("statistic" %in% names(result))
  expect_true(is.numeric(result$p.value))
  expect_true(result$p.value >= 0 && result$p.value <= 3)
})

test_that("bimax_biclust finds biclusters in a known matrix", {
  # Create a binary matrix with a known bicluster
  m <- matrix(0, nrow = 5, ncol = 5)
  m[1:3, 1:3] <- 1  # 3x3 all-ones block
  rownames(m) <- 1:5
  colnames(m) <- 1:5

  result <- BicliqueRT:::bimax_biclust(m, minr = 2, minc = 2, number = 5)
  expect_true(result$Number >= 1)
  expect_true(is.matrix(result$RowxNumber))
  expect_true(is.matrix(result$NumberxCol))
})

test_that("bimax_extract returns correct submatrices", {
  m <- matrix(0, nrow = 4, ncol = 4)
  m[1:2, 1:2] <- 1
  rownames(m) <- 1:4
  colnames(m) <- 1:4

  bic <- BicliqueRT:::bimax_biclust(m, minr = 2, minc = 2, number = 1)
  extracted <- BicliqueRT:::bimax_extract(m, bic)

  if (bic$Number > 0) {
    expect_true(is.matrix(extracted$Bicluster1))
    expect_true(all(extracted$Bicluster1 == 1))
  }
})

test_that("bimax_indices returns row and column indices", {
  m <- matrix(0, nrow = 4, ncol = 4)
  m[1:3, 1:3] <- 1
  rownames(m) <- 1:4
  colnames(m) <- 1:4

  bic <- BicliqueRT:::bimax_biclust(m, minr = 2, minc = 2, number = 1)
  indices <- BicliqueRT:::bimax_indices(bic)

  if (bic$Number > 0) {
    expect_true(is.list(indices$Bicluster1))
    expect_true("Rows" %in% names(indices$Bicluster1))
    expect_true("Cols" %in% names(indices$Bicluster1))
  }
})

test_that("greedy decomposition finds a biclique", {
  set.seed(10)
  # Create a NEgraph with some zeros so the greedy algo can partition

  NEgraph <- matrix(sample(c(0, 1), 10 * 20, replace = TRUE, prob = c(0.3, 0.7)),
                    nrow = 10, ncol = 20)
  rownames(NEgraph) <- 1:10
  colnames(NEgraph) <- 1:20

  result <- out_greedy_decom(NEgraph, num_ass = 3)
  if (!is.null(result$clique)) {
    expect_true(is.matrix(result$clique))
    expect_true(ncol(result$clique) >= 3)
  }
})

test_that("biclique.decompose works with greedy method", {
  set.seed(42)
  num_units <- 20
  Z <- rbinom(num_units, 1, 0.5)

  hypothesis <- list(
    design_fn = function() rbinom(num_units, 1, 0.5),
    exposure_i = function(z, i) z[i],
    null_equiv = function(e1, e2) identical(e1, e2)
  )

  controls <- list(method = "greedy", mina = 5, num_randomizations = 50)
  result <- suppressMessages(biclique.decompose(Z, hypothesis, controls))

  expect_true(is.list(result))
  expect_true("MNE" %in% names(result))
  expect_true(length(result$MNE) > 0)
})

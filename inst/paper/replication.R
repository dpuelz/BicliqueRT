##############################################################################
## Replication script for:
## "BicliqueRT: Randomization Tests for Causal Effects Under General
##  Interference in R"
##
## This script reproduces all examples from the JSS manuscript.
##############################################################################

library("BicliqueRT")

cat("=== BicliqueRT Replication Script ===\n\n")

## ---- Section 4.1: Spatial interference --------------------------------------

cat("--- Section 4.1: Spatial interference ---\n")

set.seed(42)
num_units <- 50

# Generate a spatial network (3 Gaussian clusters)
net <- out_example_network(num_units)
D <- net$D

# Adjacency matrix: connect units within distance 0.15
G <- (D < 0.15 & D > 0) * 1

# Treatment and outcomes
Z <- rbinom(num_units, 1, 0.3)
Y <- rnorm(num_units) + 0.5 * Z

# Hypothesis: no interference (exposure = own treatment)
hypothesis <- list(
  design_fn  = function() rbinom(num_units, 1, 0.3),
  exposure_i = function(z, i) z[i],
  null_equiv = function(e1, e2) identical(e1, e2)
)

controls <- list(method = "greedy", mina = 10, num_randomizations = 500)
decomp <- biclique.decompose(Z, hypothesis, controls)

tstat <- gen_tstat(G, type = "htn")
result <- clique_test(Y, Z, tstat, decomp)

cat("Spatial interference test p-value:", result$p.value, "\n\n")


## ---- Section 4.2: Clustered interference ------------------------------------

cat("--- Section 4.2: Clustered interference ---\n")

set.seed(123)
N <- 90; K <- 30; numrand <- 200

Zprime <- out_Zprime(N, K, equal = TRUE, numrand = numrand)
Z_a <- (Zprime == 1) * 1   # spillover
Z_b <- (Zprime == 0) * 1   # pure control
Z_treat <- (Zprime == 2) * 1

Zobs_id <- 1
bf <- out_bassefeller(N, K, Zobs = Zprime[, Zobs_id], tau_main = 0)

result2 <- clique_test_contrast(
  Y = bf$Yobs, Z = Z_treat, Z_a = Z_a, Z_b = Z_b,
  Zobs_id = Zobs_id, decom = "greedy", mina = 10
)

cat("Clustered interference test decision:", result2$decision, "\n\n")


## ---- Section 4.3: Network interference --------------------------------------

cat("--- Section 4.3: Network interference ---\n")

set.seed(99)
num_units <- 40

# Erdos-Renyi graph
G2 <- matrix(0, num_units, num_units)
for (i in 1:(num_units - 1)) {
  for (j in (i + 1):num_units) {
    if (runif(1) < 0.1) G2[i, j] <- G2[j, i] <- 1
  }
}

Z2 <- rbinom(num_units, 1, 0.3)
Y2 <- rnorm(num_units) + 0.3 * (G2 %*% Z2)

hypothesis2 <- list(
  design_fn  = function() rbinom(num_units, 1, 0.3),
  exposure_i = function(z, i) c(z[i], sum(G2[i, ] * z)),
  null_equiv = function(e1, e2) identical(e1, e2)
)

controls2 <- list(method = "greedy", mina = 5, num_randomizations = 300)
decomp2 <- biclique.decompose(Z2, hypothesis2, controls2)

tstat2 <- gen_tstat(G2, type = "score")
result3 <- clique_test(Y2, Z2, tstat2, decomp2)

cat("Network interference test p-value:", result3$p.value, "\n\n")


## ---- Session info -----------------------------------------------------------

cat("--- Session info ---\n")
sessionInfo()

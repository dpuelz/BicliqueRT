# BicliqueRT
R package for randomization tests of causal effects under general interference.  The full package is under construction, but the main functions are up and running!  

This is an implementation of the clique-based randomization test developed in the paper [A Graph-Theoretic Approach to Randomization Tests of Causal Effects Under General Interference](https://arxiv.org/pdf/1910.10862.pdf). If used, please cite the paper.

## Installation
To use, first install devtools:
```R
install.packages("devtools")
```
and then, install the package:
```R
library(devtools)
install_github("dpuelz/BicliqueRT")
```

## Example: Spatial Interference
The following simulation example illustrates spatial inference on a small synthetic network with 500 nodes:

```R
# generated network - 3 clusters of 2D Gaussians
# loads in the 500x500 matrix Dmat
# Dmat just encodes all pairwise Euclidean distances between network nodes, and
# this is used to define the spillover hypothesis below.

library(BicliqueRT)
set.seed(1)
thenetwork = out_example_network(500)
D = thenetwork$D

# simulation parameters
num_randomizations = 5000
radius = 0.01

# First, construct Z, Z_a, Z_b.
# Here, a unit receives exposure a if it is untreated, but within radius of a 
# treated unit; it receives exposure b if it is untreated and at least radius 
# distance away from all treated units.
# Experimental design is Bernoulli with prob=0.2.
# a_threshold is a scalar denoting the threshold that triggers an exposure to a.  
# If exposure a is simply binary, i.e. whether or not unit j is exposed to a, then 
# this value should be set to 1.
# b_threshold is a scalar denoting the threshold that triggers an exposure to b.  If exposure b
# is simply binary, i.e. whether or not unit j is exposed to b, then this value should be set to 1.
Z = out_Z(pi=rep(0.2,dim(D)[1]),num_randomizations)
D_a = D_b = sparsify((D<radius))
a_threshold = b_threshold = 1

Z_a = Z_b = D_a%*%Z
Z_a = sparsify((Z_a>=a_threshold))
Z_b = sparsify((Z_b<b_threshold))

# simulating an outcome vector
Y = rnorm(dim(Z)[1])

# run the test using Bimax to decompose the null-exposure graph
# tau=0.2 is the tau in the null: Y_i(b) = Y_i(a) + tau for all i.
CRT = clique_test(Y, Z, Z_a, Z_b, Zobs_id=1, tau=0.2, decom='bimax', minr=15, minc=15)

# alternatively, we can use a greedy algorithm to do decomposition by specifying decom
CRT = clique_test(Y, Z, Z_a, Z_b, Zobs_id=1, tau=0.2, decom='greedy', minass=15)
```
Sometimes we want to replace the outcome vector Y with an adjusted version. We can pass in Xadj and specifying adj_Y=TRUE. Currently we only support adjusing Y by taking the residuals of a linear regression on Xadj, but users can also pre-adjust it before using the clique test function.
```R
Xadj = matrix(rnorm(dim(Z)[1]*4), ncol=4)
CRT = clique_test(Y, Z, Z_a, Z_b, Zobs_id=1, Xadj=Xadj, tau=0.2, decom='bimax', minr=15, minc=15, adj_Y=TRUE)
```
To get the CI, we can pass in ret_ci=TRUE, and it's not necessary to specify tau in this case:
```R
CRT = clique_test(Y, Z, Z_a, Z_b, Zobs_id=1, decom='bimax', ret_ci=TRUE, ci_method='grid', minr=15, minc=15)
```
By default, we use the "grid" method to calculate CI, we provide another method "bisection". 
We can also do parallelization for the grid method by specifying Cluster beforehand:

```R
library(doParallel)
numcores = detectCores()
clst = makeCluster(numcores[1]-1) 
registerDoParallel(clst)
CRT = clique_test(Y, Z, Z_a, Z_b, Zobs_id=1, decom='bimax', ret_ci=TRUE, ci_method='grid', minr=15, minc=15)
```

## Example: Clustered Interference
The following simulation example illustrates clustered inference with 2000 individuals equally divided into 500 clusters:

```R
library(BicliqueRT)
set.seed(1)
N = 2000 # total number of individuals
K = 500  # total number of households, i.e., number of clusters
Zobs_id = 1

# Generate household-individual structure and experiment design.
# Each column of Zprime_mat specifies an assignment, and each row represents an individual.
# Entries of Zprime_mat is either 0, 1, or 2, indicating individual's exposure.
# Here, an individual has exposure 0 if it's within untreated cluster;
# it has exposure 1 if it's untreated but someone else in the same cluster is treated;
# it has exposure 2 if it's treated.
Zprime_mat = out_Zprime(N, K, numrand=1000)
Z = Zprime_mat==2    # we convert Z to be a binary matrix indicating whether individual is treated (T) or not (F)
Z_a = Zprime_mat==1  # controlled individuals in treated households, "spillover"
Z_b = Zprime_mat==0  # individuals in untreated households, "controlled"

# simulate an outcome vector assuming the null is true
simdat = out_bassefeller(N, K, Zprime_mat[, Zobs_id],tau_main = 0.4)
Yobs = simdat$Yobs

# run the test using Bimax
CRT = clique_test(Yobs, Z, Z_a, Z_b, Zobs_id, tau=0, decom='bimax', minr=20, minc=20)

# again, we can use the greedy algorithm as follows:
CRT = clique_test(Yobs, Z, Z_a, Z_b, Zobs_id, tau=0, decom='greedy', minass=20)

```

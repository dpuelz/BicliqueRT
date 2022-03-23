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
Sometimes we want to replace the outcome vector Y with an adjusted version. We can pass in Xadj and specifying adj_Y=TRUE. Currently we only support adjusting Y by taking the residuals of a linear regression on Xadj, but users can also pre-adjust it before using the clique test function.
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
clst = makeCluster(numcores-1) # try not to use all cores in order not to impose great burden on the computer
registerDoParallel(clst)
CRT = clique_test(Y, Z, Z_a, Z_b, Zobs_id=1, decom='bimax', ret_ci=TRUE, ci_method='grid', minr=15, minc=15)
stopCluster(clst)
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
simdat = out_bassefeller(N, K, Zprime_mat[, Zobs_id], tau_main = 0.4)
Yobs = simdat$Yobs
```
This simulate the cluster structure data such that individuals in spillover groups (house i treated and self not treated) are on average 0.4 higher than individuals in pure controlled groups (house i not treated and self not treated).
```R
# run the test under the null that tau = 0
CRT = clique_test(Yobs, Z, Z_a, Z_b, Zobs_id, tau=0, decom='bimax', minr=25, minc=25)
CRT$decision
# [1] 0
```
we get literally 0 p-value, which strongly rejects the null that there is no spillover effect. 

If we instead specify the null to be that the spillover effect is exactly 0.4 in the test:
```R
# run the test under the null that tau = 0.4
CRT1 = clique_test(Yobs, Z, Z_a, Z_b, Zobs_id, tau=0.4, decom='bimax', minr=25, minc=25)
CRT1$decision
# [1] 0.9166667
```
or if we change the simulated data such that there is actually on average no difference between the two groups:
```R
simdat2 = out_bassefeller(N, K, Zprime_mat[, Zobs_id], tau_main = 0)
Yobs2 = simdat2$Yobs
# run the test under the null that tau = 0
CRT2 = clique_test(Yobs2, Z, Z_a, Z_b, Zobs_id, tau=0, decom='bimax', minr=25, minc=25)
CRT2$decision
# [1] 0.5
```
we can see that in both cases, we cannot reject the null, as desired.

## Example: Diagonistic graph
A common clustered interference problems has the following structure of the experiment design:
- Each group has a probability of `p_group` to be selected as the treated group
- Individuals in the treated groups have a probability of `p_indi_t` to be treated, and individuals in the untreated groups have a probability of `p_indi_nt` to be treated

And we would like to test whether untreated individuals in untreated groups have similar outcomes to untreated individuals in treated groups. For example, see [Breza et al., 2021](https://www.nature.com/articles/s41591-021-01487-3).

Usually the null exposure graph of such experiment design is very "spare", in the sense that there are too few exposures to construct a large enough biclique. The `clique_diagnostic` function can serve as a diagnostic for how large we could expect the biclique could be. The input of the function is a (# of individuals) * 2 matrix or data frame whose first column indicates group, and second column indicates individual ID. Given the structure of the group-individual, for different number of individuals selected at random, the function returns the approximate number of focal assignments of the largest possible biclique of the null exposure graph that takes the selected units as focal units. 

Below is an example illustrating the function. The group-individual structure used is the Thanksgiving campaign of [Breza et al., 2021](https://www.nature.com/articles/s41591-021-01487-3).

```R
data("sample_structure")
sample_structure[1:3,] # user_loc is county code, zip is zip code
#  user_loc   zip
#1     4001 85938
#2     4001 85925
#3     4001 86505

set.seed(1)
# we try for 10000 randomizations instead of 1000 in the original article.
clique_diag = clique_diagnostic(struc = sample_structure, p_group = 0.5, p_indi_t = 0.75, p_indi_nt = 0.25, N = 10000) # takes about 10 min to complete
t(clique_diag)
#        [,1]      [,2]      [,3]     [,4]     [,5]     [,6]     [,7]     [,8]     [,9]    [,10]
#[1,] 5000.82 2508.7250 1255.5637 626.1600 313.3240 157.8799 79.34066 39.60684 20.07522 9.919738
#[2,]    0.00  934.8375  699.8063 423.9984 237.0436 128.1890 68.26418 35.54239 18.40697 9.329229
#        [,11]    [,12]    [,13]     [,14]     [,15]     [,16]      [,17]      [,18]      [,19]
#[1,] 5.107785 2.528976 1.271910 0.6466321 0.3259316 0.1616551 0.08285015 0.04274067 0.02168420
#[2,] 4.875584 2.451729 1.243085 0.6330025 0.3175506 0.1586706 0.08194811 0.04203073 0.02151507
#          [,20]       [,21]       [,22]      [,23]        [,24]       [,25]
#[1,] 0.01067952 0.005542267 0.002781877 0.00139669 0.0007813702 0.000354051
#[2,] 0.01055267 0.005531696 0.002781877 0.00139669 0.0007813702 0.000354051
```
we can also plot it to see how fast it declines when we increase the number of units selected as focal units.
```R
plot(x = 2:25, y = clique_diag[-1,1], 'l', lty = 1, col = 'red',
     xlab = "num focal units", ylab = "num of focal assin ", ylim = c(0,30)) # total number of focal assignments
lines(x = 2:25, y = clique_diag[-1,2], 'l', lty = 1, col = 'green') # number of focal assignments that have variations across focal units
```
So if we want the final biclique decomposed from the null exposure graph that contains `Zobs` to have at least 12 focal units, and also each of the biclique's focal assignments does not contain only one type of exposure (so that we can do randomization test on the biclique), it is likely to be very small because on average it contains only 2.45 focal assignments as indicated above. If in the `Bimax` algorithm we set `minr = 12` and, say, `minc = 15`, it would take quite a long time to decompose the null exposure graph. What's worse is that it may never find such a biclique that contains `Zobs`!

## Example: Extent of Interference
We demonstrate how to test the "extent of interference" type null hypotheses as in example 3 in the paper. The test allows the individual exposure to be multi-dimensional, and test whether for all individuals, the potential outcome is the same under different treatment assignments that give the same exposure for an individual (eq. 4 in the paper). We illustrate using the example 3 in the paper.

We firstly generate the network
```R
set.seed(1)
N = 30
D = matrix(sample(c(1:3, Inf), N^2, prob = c(.4, .3, .2, .1), replace = T), N, N)
D[lower.tri(D)] = t(D)[lower.tri(D)]
diag(D) = 0
```
`D` is a `N` by `N` symmetric matrix where each measures the distance between units `i` and `j` in the network. The distance takes five values 0,1,2,3,Inf where `D[i,i]=0`, and `D[i,j]=Inf` if `i` and `j` are not connected. Smaller the value, closer unit `i` and `j` are.

The treatment assignment mechanism is an individual Bernoulli trial with treated probability being 0.2. 
```R
num_randomizations = 1000
Z = out_Z(pi=rep(0.2, dim(D)[1]), num_randomizations)
```
We let `k=1` that assumes individuals' potential outcomes may depend only on treatments of units up to 1 hops away in the network, but no further.




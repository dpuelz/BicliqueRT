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



## General Work Flow of the Test

Suppose we are going to test the null: Yi(z)=Yi(z') for all i and all z,z' such that fi(z)~fi(z'), where fi(z) is the exposure of unit i under treatment assignment z, and fi(z)~fi(z') expresses the equivalence of two exposures.

The test consists of two parts as described in the paper: do biclique decomposition on the null exposure graph to find a conditional clique, and do randomization test on the conditional clique.

### 1. Biclique decomposition

The function in this step is `biclique.decompose(Z, hypothesis, controls=list(method="greedy", mina=10, num_randomizations=2000), stop_Zobs=F)`. `Z` is the observed treatment assignment vector of length N where N is the number of units. `hypothesis` is a list that contains three functions specifying the experiment design and null hypothesis: `design_fn`, `exposure_i` and `null_equiv`.

* `design_fn` specifies the experiment design. It should be a function that returns a realization of treatment assignment for the whole sample. It can depends on other global variables such as the number of units. For example, if the design is that each unit has equal probability of 0.2 to receive the treatment independently, we can write `design_fn = function() { rbinom(num_units, 1, prob=0.2) }`. Note that the return of `design_fn()` should also be a vector of length N as `Z` is.
* `exposure_i` should be a function that returns exposure fi(z) of unit `i` under treatment `z` for the whole sample. The inputs of the function are an index `i` and a vector `z`. For example, if the exposure of `i` under `z` is the treatment it receives, then we can write `exposure_i = function(z, i) { z[i] }`. See examples below for more instances of `exposure_i`.
* `null_equiv` expresses the equivalence relationship between two exposures. It should be a function that takes two inputs from `exposure_i `and determines whether they are equivalent according to the null hypothesis. See examples below for more instances of `null_equiv`.

The list `controls` specifies parameters for the biclique decomposition algorithm. `method` could be either `"bimax"` or `"greedy"` that uses the two biclique decomposition algorithm. If `"bimax"` is used, `minr` and `minc ` should be supplied that specify the minimum number of units and assignments in the bicliques found by the algorithm. If `"greedy"` is used, `mina ` should be supplied. `num_randomizations` specifies the number of randomizations to perform. Larger the number, more computation time is needed but we may get a larger conditional clique to do randomization which gives the test more power.

`stop_Zobs` is either `TRUE` or `FALSE` that determines whether the biclique decomposition should end when we find a biclique that contains the observed treatment assignment vector. Setting it to `TRUE` can speed up the decomposition a bit but may not get the whole biclique decomposition of the null exposure graph.

The `biclique.decompose` function will return a (partial) biclique decomposition of the null exposure graph where `Z` is in one of them.

### 2. Randomization Test

The output of the `biclique.decompose`, written as `biclique_decom`, should be passed to `clique_test(Y, Z, teststat, biclique_decom, alpha=0.05)` to do the randomization test. Here `Y` is the observed outcome vector of length N, `Z` is the observed treatment assignment vector of same length. `alpha` specifies the significance level.

`teststat` is a function that specifies the test statistic used in the conditional clique. The function should contain at least (with order) `y, z, focal_unit_indicator` as inputs, where `y` is the outcome vector, `z` is the treatment vector and `focal_unit_indicator` is a 0-1 vector indicating whether a unit is focal (=1) or not (=0). All three inputs should have length equal to number of units and have the same ordering. Other global variables can be used in the function, such as the ones about a network of interference.

We provide several default test statistics for no-interference null in [Athey et al. 2018 Exact p-Values for Network Interference](https://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1241178) that can be generated using function `gen_tstat(G, type)`, where `G` is the N by N adjacency matrix of a network, `type` could be one of `"elc","score","htn"`. See section 5 in [Athey et al. 2018](https://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1241178) and the documentation of `gen_tstat` for a detailed description of the these test statistics.



## Example: Spatial Interference

The following simulation example illustrates spatial inference on a small synthetic network with 500 nodes. We are testing the null that potential outcomes are the same for all units no matter (it is untreated but within a radius of a treated unit), or (it is untreated but not within a radius of any treated unit).

```R
# generated network - 3 clusters of 2D Gaussians
# loads in the 500x500 matrix Dmat
# Dmat just encodes all pairwise Euclidean distances between network nodes, and
# this is used to define the spillover hypothesis below.

library(BicliqueRT)
set.seed(1)

N = 500 # number of units
thenetwork = out_example_network(N)
D = thenetwork$D

# simulating an outcome vector and a treatment realization
Y = rnorm(N)
Z = rbinom(N, 1, prob=0.2)

# simulation parameters
num_randomizations = 5000
radius = 0.02

# To use the package:
#   1. The design function: here the experimental design is Bernoulli with prob=0.2
design_fn = function() { rbinom(N, 1, prob=0.2) }

# 	2. The exposure function: exposure for each unit is (w_i, z_i) where
# 			w_i = 1{\sum_{j\neq i} g_{ij}^r z_j > 0 } and g_{ij}^r = 1{d(i,j)<r}
Gr = (D<radius) * 1; diag(Gr) = 0
exposure_i = function(z, i) { c(as.numeric(sum(Gr[i,]*z) > 0), z[i]) }

# 	3. The null
null_hypothesis = list(c(0,0), c(1,0))
null_equiv = function(exposure_z1, exposure_z2) {
  (list(exposure_z1) %in% null_hypothesis) & (list(exposure_z2) %in% null_hypothesis)
}

# Then we can decompose the null exposure graph:
H0 = list(design_fn=design_fn, exposure_i=exposure_i, null_equiv=null_equiv)
bd = biclique.decompose(Z, H0, controls= list(method="greedy", mina=20, num_randomizations = 2e3))
m = bd$MNE # this gives the biclique decomposition

# To do randomization test, firstly generate a test statistic. Here we use the contrast of mean between units with exposure (0,0) and exposure (1,0)
Tstat = function(y, z, is_focal) {
  exposures = rep(0, N)
  for (unit in 1:N) {
    exposures[unit] = exposure_i(z, unit)[1]
  }
  stopifnot("all focals have same exposures" = (length(unique(exposures[is_focal]))>1) )
  mean(y[is_focal & (exposures == 1)]) - mean(y[is_focal & (exposures == 0)])
}

# Then run the test
testout = clique_test(Y, Z, Tstat, bd)
```
The output `testout` is a list that contains p-value, the test statistic, the randomization distribution of the test statistic, test method and the conditional clique.

```R
names(testout)
# [1] "p.value"            "statistic"          "statistic.dist"     "method"             "conditional.clique"
testout$p.value # p-value of the test
```



## Example: Clustered Interference

The following simulation example illustrates clustered inference with 2000 units equally divided into 500 clusters. We assign at random one unit in half of the clusters to be treated, and we are testing that potential outcomes are the same for all units no matter (it is untreated in a cluster without any treated unit), or (it is untreated in a cluster with a treated unit). We follow the same procedure:

```R
library(BicliqueRT)
set.seed(1)
N = 2000 # total number of units
K = 500  # total number of households, i.e., number of clusters
housestruct = out_house_structure(N, K, T)

# The design function:
design_fn = function(){
  treatment_list = out_treat_household(housestruct, K1 = K/2) # one unit from half of the households would be treated.
  treatment = out_Z_household(N, K, treatment_list, housestruct)
  return(treatment[,'treat'])
}

# Generate a treatment realization and outcome
Z = design_fn() # randomly get one realization
Y = out_bassefeller(N, K, Z, tau_main = 0.4, housestruct = housestruct)$Yobs 
# here we assume that potential outcomes are 0.4 higher if an untreated unit is in a cluster
# with a treated unit compared to in a cluster without any treated unit, 
# i.e., a spillover effect of 0.4 is assumed

# The exposure function: exposure for each unit i is z_i + \sum_{j \in [i]} z_j where [i] represents the cluster i is in.
exposure_i = function(z, i) {
  # find the household that i is in
  house_ind = cumsum(housestruct)
  which_house_i = which.min(house_ind < i)
  # find lower and upper index of [i] in z
  if (which_house_i == 1) {lower_ind = 1} else {lower_ind = house_ind[which_house_i-1] + 1}
  upper_ind = house_ind[which_house_i]
  # calculate exposure
  exposure_z = z[i] + sum(z[lower_ind:upper_ind])
  exposure_z
}

# The null
null_equiv = function(exposure_z1, exposure_z2) {
  ((exposure_z1 == 1) | (exposure_z1 == 0)) & ((exposure_z2 == 1) | (exposure_z2 == 0))
}

# Do biclique decomposition on the null exposure graph
H0 = list(design_fn=design_fn, exposure_i=exposure_i, null_equiv=null_equiv)
bd = biclique.decompose(Z, H0, controls= list(method="greedy", mina=30, num_randomizations = 5e3))

# Define a test statistic, here we use the contrast of mean between units with exposure 1 (untreated in treated cluster) and exposure 0 (untreated in untreated cluster)
Tstat = function(y, z, is_focal) {
  exposures = rep(0, N)
  for (unit in 1:N) {
    exposures[unit] = exposure_i(z, unit)[1]
  }
  stopifnot("all focals have same exposures" = (length(unique(exposures[is_focal]))>1) )
  mean(y[is_focal & (exposures == 1)]) - mean(y[is_focal & (exposures == 0)])
}

# Then run the test
testout = clique_test(Y, Z, Tstat, bd)
testout$p.value # p-value of the test
```



## Example: Diagnostic graph

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
So if we want the final biclique decomposed from the null exposure graph that contains the observed `Z` to have at least 12 focal units, and also each of the biclique's focal assignments does not contain only one type of exposure (so that we can do randomization test on the biclique), it is likely to be very small because on average it contains only 2.45 focal assignments as indicated above. If in the `Bimax` algorithm we set `minr = 12` and, say, `minc = 15`, it would take quite a long time to decompose the null exposure graph. What's worse is that it may never find such a biclique that contains the observed `Z`!



## Example: Extent of interference
Consider a setting where units are linked in a network. The null we are testing is whether for all units, the potential outcome depends only on treatments of units that are k distant from the unit itself. We illustrate using the case of k=0, which is the null of no interference that a unit's potential outcome depends only on its treatment. 

We firstly generate the network and potential outcomes. We use the Erdos-Renyi model where every possible edge is created with the same constant probability 0.1. We assign exactly half of the units to be treated.
```R
library(BicliqueRT)
library(igraph)
set.seed(1)

N = 100
g = erdos.renyi.game(N, 0.1)

Z = sample(c(rep(1, N/2), rep(0, N/2))) # treatment
A = get.adjacency(g)
G = as.matrix(A); diag(G) = 0

W = as.numeric(G %*% Z)
Y = 0.1 + 1*Z  + 5*W+ rnorm(N) # here we have a positive spillover effect from immediate neighbours
```
`G` is a `N` by `N` symmetric matrix where `G[i,i]=0`, and `G[i,j]=1` if there's an edge between units `i` and `j` in the network. We then do the same things to decompose the null exposure graph.

~~~R
# The design function
design_fn = function() sample(c(rep(1, N/2), rep(0, N/2)))

# The exposure function
exposure_i = function(z, i){
  stopifnot(length(z)==N)
  z[i]
}

# The null
null_equiv = function(e1, e2){
  identical(e1, e2)
}

# decompose
H0 = list(design_fn=design_fn, exposure_i=exposure_i, null_equiv=null_equiv)
bd = biclique.decompose(Z, H0, controls= list(method="greedy", mina=50, num_randomizations = 2e3))
~~~

We can use the `gen_tstat` function to generate a test statistic. Here we use the edge-level-contrast statistic. Then we can run the test.

```R
Tstat = gen_tstat(G, "elc")
testout = clique_test(Y,Z, Tstat, bd)
testout$p.value # p-value of the test
```

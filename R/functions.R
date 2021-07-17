#' The main randomization test function.
#'
#' @param Y The observed outcome vector.
#' @param Z A binary matrix of dimension (number of units x number of randomizations, i.e. assignments.) storing the assignment vectors. Please see example.
#' @param Z_a A binary matrix with dimension (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{a} under assignment j. Please see example.
#' @param Z_b A binary matrix with (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{b} under assignment j. Please see example.
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.
#' @param minr The minimum number of focal units included in the cliques (algorithm runtime is sensitive to this value).
#' @param minc The minimum number of focal assignment included in the cliques (algorithm runtime is sensitive to this value).
#' @param ... Other stuff ...
#'
#' @return A list of items summarizing the randomization test.
#' @examples
#' # Spatial interference
#' # generated network - 3 clusters of 2D Gaussians
#' # loads in the 500x500 matrix Dmat
#' # Dmat just encodes all pairwise Euclidean distances between network nodes, and
#' # this is used to define the spillover hypothesis below.
#' library(CliqueRT)
#' set.seed(1)
#' thenetwork = out_example_network(500)
#' D = thenetwork$D
#'
#' # simulation parameters
#' num_randomizations = 5000
#' radius = 0.01
#'
#' # First, construct Z, Z_a, Z_b.
#' # Here, exposure a is an untreated within radius of a treated unit, and exposure b
#' # is an untreated unit at least radius distance away from all treated units.
#' # Experimental design is Bernoulli with prob=0.2.
#' # a_threshold is a scalar denoting the threshold that triggers an exposure to a.  If exposure a
#' # is simply binary, i.e. whether or not unit j is exposed to a, then this value should be set to 1.
#' # b_threshold is a scalar denoting the threshold that triggers an exposure to b.  If exposure b
#' # is simply binary, i.e. whether or not unit j is exposed to b, then this value should be set to 1.
#' Z = out_Z(pi=rep(0.2,dim(D)[1]),num_randomizations)
#' D_a = D_b = sparsify((D<radius))
#' a_threshold = b_threshold = 1
#'
#' Z_a = Z_b = D_a%*%Z
#' Z_a = sparsify((Z_a>=a_threshold))
#' Z_b = sparsify((Z_b<b_threshold))
#'
#' # simulating an outcome vector
#' Y_a = rnorm(dim(Z)[1])
#' Y_b = Y_a + 0.2
#' Y = out_Yobs(Z_a[,1],Z_b[,1],Y_a,Y_b)
#'
#' # run the test
#' CRT = clique_test(Y,Z,Z_a,Z_b,Zobs_id=1,minr=15,minc=15)
#'
#' # Clustered interference
#' # simulation parameters
#' N = 2000
#' K = 500
#' Zobs_id = 1
#'
#' # generate clustered structure
#' library(Matrix)
#' library(biclust)
#' set.seed(1)
#' Zprime_mat = out_Zprime(N, K, numrand=1000)
#' Z = Zprime_mat==2
#' Z_a = Zprime_mat==1
#' Z_b = Zprime_mat==0
#'
#' # simulate an outcome vector
#' simdat = out_bassefeller(N, K, Zprime_mat[, Zobs_id],tau_main = 0.4)
#' Yobs = simdat$Yobs
#'
#' # run the test
#' CRT = clique_test(Yobs, Z, Z_a, Z_b, Zobs_id, minr = 20, minc = 20)
#' @export
clique_test = function(Y,Z,Z_a,Z_b,Zobs_id,minr,minc,...){

  # setting default values if they do not exist
  if(!exists("exclude_treated")){ exclude_treated=TRUE }
  if(!exists("stop_at_Zobs")){ stop_at_Zobs=FALSE }
  if(!exists("ret_pval")){ ret_pval=TRUE }

  # make the null-exposure graph
  cat("construct the null-exposure graph ... \n")
  NEgraph = out_NEgraph(Z_a,Z_b,Z,exclude_treated)

  # decompose the null-exposure graph
  cat("decompose the null-exposure graph ... \n")
  decomp = out_clique_decomposition(NEgraph,Zobs_id,minr,minc,stop_at_Zobs)
  if(stop_at_Zobs){ conditional_clique = decomp }
  if(!stop_at_Zobs){ conditional_clique = out_clique(Zobs_id,decomp) }

  focal_units = as.numeric(rownames(conditional_clique)) # a list of row names
  focal_assignments = as.numeric(colnames(conditional_clique)) # a list of column names

  # restrict Y to clique
  Y.clique = Y[focal_units]

  # run the test
  cat("\n")
  cat("run the clique-based randomization test ... \n")
  Zobs_cliqid = which(Zobs_id==focal_assignments)
  tobs = ate(conditional_clique[, Zobs_cliqid], Y.clique)
  tvals = apply(as.matrix(conditional_clique[, -Zobs_cliqid]), 2, function(z) ate(z, Y.clique))
  rtest_out = list(tobs=tobs,tvals=tvals)
  decision = out_pval(rtest_out,ret_pval,alpha=0.05)

  # organize into return list
  retlist = list(decision=decision,ret_pval=ret_pval,tobs=tobs,tvals=tvals,focal_units=focal_units,focal_assignments=focal_assignments,NEgraph=NEgraph)
  retlist
}


#' One-sided testing
#'
#' Decides to reject or not based on observed test statistic value \code{tobs} and randomization values \code{tvals}.
#'
#' @param tobs The observed value of the test statistic (scalar).
#' @param tvals Vector of randomization values of the test statistic (to compare with \code{tobs}).
#' @param alpha Desired level of the test (between 0 to 1).
#' @param tol Used to check whether \code{tobs} is equal to the 1-\code{alpha} quantile of \code{tvals}.
#' If the observed \code{tobs} is within \code{tol} of any of \code{tvals}, they
#' will be treated as equal.
#' @details
#' The test may randomize to achieve the specified level \code{alpha}
#' when there are very few randomization values. Returns 1 if the test rejects, 0 otherwise.
#'
#' Note that if the \eqn{1-\alpha} percentile of \code{tvals} is the same as
#' \code{tobs} (up to \code{tol}), it will return a randomized decison.
#' @return Test decision (binary). Returns 1 if the test rejects, 0 otherwise.
#' @examples
#' tvals <- seq(0.1, 1, length.out=10)
#' tobs <- 0.85
#' > one_sided_test(tobs, tvals, 0.1)
#' [1] 0
#' @seealso Testing Statistical Hypotheses (Ch. 15, Lehman and Romano, 2006)
#' @export
one_sided_test = function(tobs, tvals, alpha, tol=1e-14) {
  srt = sort(tvals)
  M = length(tvals)
  k = ceiling(M * (1-alpha))
  Tk = srt[k]
  if(abs(tobs - Tk) < tol) {
    # if tobs = Tk
    ax = (M * alpha - sum(tvals > Tk)) / sum(abs(tvals - Tk) < tol)
    return(1*(runif(1) <= ax)) ## randomize decision.
  }

  return(1*(tobs > Tk))
}

#' Two-sided testing
#'
#' Decides to reject or not based on observed test statistic value \code{tobs}
#' and randomization values \code{tvals}. The test may randomize to achieve the specified level \code{alpha}
#' when there are very few randomization values.
#'
#' @param tobs The observed value of the test statistic (scalar).
#' @param tvals Vector of randomization values of the test statistic (to compare with \code{tobs}).
#' @param alpha Desired level of the test (between 0 to 1).
#' @param tol Used to check whether \code{tobs} is equal to the 1-\code{alpha}/2
#' or \code{alpha}/2 quantile of \code{tvals}.
#' If the observed \code{tobs} is within \code{tol} of any of \code{tvals}, they
#' will be treated as equal.
#'
#' @examples
#' tvals <- seq(0, 1, length.out=1001)
#' tobs <- 0.95 + 1e-13
#' > two_sided_test(tobs, tvals, 0.1)
#' [1] 1
#'
#' @return Test decision (binary). Returns 1 if the test rejects, 0 otherwise.
#' @seealso Testing Statistical Hypotheses (Ch. 15, Lehman and Romano, 2006)
#' @export
two_sided_test = function(tobs, tvals, alpha, tol=1e-14) {
  m1 = one_sided_test(tobs, tvals, alpha=alpha/2, tol=tol)
  m2 = one_sided_test(-tobs, -tvals, alpha=alpha/2, tol=tol) # only one can be 1.
  return(m1 + m2)
}

#' Calculates p-value or test decision
#'
#' Depending on \code{ret_pval} this function returns either a p-value for the test or the binary decision.
#'
#' @param rtest_out A \code{List} with elements \code{tobs}, \code{tvals} (see \link{one_sided_test} for details.)
#' @param ret_pval A \code{Boolean} indicating whether to return a p-value (TRUE) or not.
#' @param alpha Desired test level (from 0 to 1).
#' @param tol tolerance level for equality between \code{tobs} and \code{tvals}.
#' @return Binary decision if \code{ret_pval} is TRUE, or the p-value otherwise.
#' @details Returns 1 if the test rejects, 0 otherwise. Note that the test is a two-sided one.
#' @export
out_pval = function(rtest_out, ret_pval, alpha, tol = 1e-14) {
  tobs = rtest_out$tobs
  tvals = c(rtest_out$tvals)

  n_all = length(tvals)
  n_higher = sum(tvals > (tobs + tol))
  n_lower = sum(tvals < (tobs - tol))
  n_equal = n_all - n_lower - n_higher

  p1 = (n_equal + n_higher) / n_all  # P(T >= Tobs)
  p2 = (n_equal + n_lower) / n_all  # P(T <= Tobs)

  pval = 2*min(p1, p2) ## since it's a two-sided test, the 2*min gives a less conservative p-value.
  if(ret_pval) return(pval)
  return(two_sided_test(tobs, tvals, alpha = alpha))  # this test is less conservative.
}


#' Computing observed outcome
#'
#' Computes the observed outcome vector, given potential outcomes on the exposures
#'
#' @param z_a A binary vector indicating which units are exposed to \code{a}.
#' @param z_b A binary vector indicating which units are exposed to \code{b}.
#' @param Y_a The potential outcome vector for all units under exposure \code{a}.
#' @param Y_b The potential outcome vector for all units under exposure \code{b}.
#'
#' @return A single observed outcome vector.
#'
#' @examples
#' Za <- c(1,1,0)
#' Zb <- c(0,0,1)
#' # Za and Zb tells us that the first and second individual are exposed to \code{a},
#' # while the third individual is exposed to \code{b}.
#' Ya <- rep(0,3)
#' Yb <- Ya + 1
#' # Ya and Yb tells us that all individuals have potential outcome 0 if exposed to \code{a}, 1 if exposed to \code{b}
#' out_Yobs(Za, Zb, Ya, Yb) # this is what we observed in this trial.
#' [1] 0 0 1
#' @export
out_Yobs = function(z_a,z_b,Y_a,Y_b){
  y = Y_a*z_a + Y_b*z_b
  y
}


#' Generating an experiment design
#'
#' Returns sparse matrix \code{Z} of dimension (number of units x number of randomizations, i.e. assignments).
#' Each column represents a randomization, with equal probability. Within each
#' randomization, individuals have a probability of \code{pi} of being treated (equals 1).
#'
#' @param pi A \code{vector} of length of number of units, specifying the
#' probability of an individual being treated in a randomization.
#' @param num_randomizations A \code{scalar} denoting the number of assignments desired.
#'
#' @return \code{Z} (sparse binary matrix). The first column is the observed assignment.
#' @export
out_Z = function(pi,num_randomizations){
  num_units = length(pi)
  Z = Matrix(rbinom(num_randomizations*num_units,1,prob=pi),nrow=num_units,ncol=num_randomizations,sparse=TRUE)
  Z
}


#' Generating Null Exposure Graph
#'
#' One of the main functions for implementing the methodology.  Outputs the null-exposure graph based on binary matrices describing exposure conditions in the null hypothesis.  The null hypothesis is represented as:
#' \eqn{H_0}: \eqn{Y_i}(\code{a}) \eqn{= Y_i}(\code{b}) for all \eqn{i},
#' and states that potential outcomes are equal for all units exposed to either \code{a} or \code{b}.
#'
#' @param Z_a A binary matrix with dimension (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{a} under assignment j.  Please see example.
#' @param Z_b A binary matrix with (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{b} under assignment j.  Please see example. \code{Z_b} should have the same dimension as \code{Z_a}.
#' @param Z A binary matrix of dimension (number of units x number of randomizations, i.e. assignments.) storing the assignment vectors.
#' @param exclude_treated A Boolean denoting whether or not treated units are considered in the hypothesis.  Default is \code{TRUE}.
#'
#' @return \code{NEgraph}, a matrix of dimension (number of units x number of randomizations, i.e. assignments.).  Row i, column j of the matrix either equals -1 (unit i is exposed to \code{a} under assignment j), 1 (unit i is exposed to \code{b} under assignment j), 0 (unit i is neither exposed to \code{a} nor \code{b} under assignment j).
#' @export
out_NEgraph = function(Z_a,Z_b,Z,exclude_treated=TRUE){
  Z_a = sparsify(Z_a)
  Z_b = sparsify(Z_b)
  Z = sparsify(Z)

  # if exclude_treated==TRUE, then component multiply Z_a, Z_b by !Z.
  if(exclude_treated){
    Z_a = Z_a*(!Z)
    Z_b = Z_b*(!Z)
  }

  # construct NEgraph from the exposures matrices
  NEgraph = Matrix(0,nrow=dim(Z_a)[1],ncol=dim(Z_a)[2],sparse=TRUE)
  NEgraph[Z_a==1] <- -1
  NEgraph[Z_b==1] <- 1

  # return
  rownames(NEgraph) = 1:nrow(NEgraph)
  colnames(NEgraph) = 1:ncol(NEgraph)
  as.matrix(NEgraph)
}


#' Finding conditioned clique
#'
#' Identifies the clique to condition on in the randomization test. Specifically, given the observed treatment
#' and a biclique decomposition, it returns a submatrix of \code{NEGraph} that represents
#' a biclique in the decomposition which contains the observed treatment.
#'
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.
#' @param decomp A list containing the clique decomposition of the null-exposure graph, given by the function \code{out_clique_decomposition}.
#'
#' @return The clique to condition on in the randomization test.
#' @export
out_clique = function(Zobs_id,decomp){
  for(ii in 1:length(decomp)){
    if( sum( Zobs_id==colnames(decomp[[ii]]) ) ){
      return(as.matrix(decomp[ii][[1]]))
    }
  }
}


#' Decomposing Null Exposure Graph NOT SURE WHY !=0
#'
#' One of the main functions for implementing the methodology.
#' Outputs a biclique decomposition of the null-exposure graph.
#' Specifically, the resulting decomposition paritions the assignment space,
#' while the unit space may overlap.
#'
#' @param NEgraph The null-exposure graph object, see \code{out_NEgraph}.
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.
#' @param minr The minimum number of focal units included in the cliques (algorithm runtime is sensitive to this value).
#' @param minc The minimum number of focal assignment included in the cliques (algorithm runtime is sensitive to this value).
#' @param stop_at_Zobs A Boolean indicating whether the decomposition algorithm should stop when the clique containing the observed assignment is found.  Default value is \code{TRUE}.
#'
#' @return If \code{stop_at_Zobs} is \code{TRUE}, a matrix representing the clique to condition upon. If \code{stop_at_Zobs} is \code{FALSE}, a list containing the clique decomposition of the null-exposure graph.
#' @export
out_clique_decomposition = function(NEgraph,Zobs_id,minr,minc,stop_at_Zobs=TRUE){

  iremove = which(rowSums(NEgraph!=0)==0)  # removes isolated units.
  if(length(iremove)!=0){ NEgraph = NEgraph[-iremove,] }
  numb = 1
  ncol = ncol(NEgraph)

  decomp = c()
  new.NEgraph = NEgraph[,1:ncol]
  numleft = ncol(new.NEgraph)
  oldnames = colnames(new.NEgraph)
  ii=1

  while(numleft>0){
    minc.new = min(minc,numleft)
    bitest = biclust(new.NEgraph, method=BCBimax(), minr, minc.new, number=numb)
    bicliqMat = bicluster(new.NEgraph,bitest)
    themat = bicliqMat$Bicluster1
    while(length(themat)==0){   ### if current minc gives no biclique decomposition, try a smaller one.
      numleft=numleft-1
      minc.new = min(minc,numleft)
      bitest = biclust(new.NEgraph!=0, method=BCBimax(), minr, minc.new, number=numb)
      bicliqMat = bicluster(new.NEgraph,bitest)
      themat = bicliqMat$Bicluster1
    }
    dropnames = colnames(themat)
    oldnames = colnames(new.NEgraph)
    drop.ind = match(dropnames,oldnames)

    if(stop_at_Zobs){
      if(sum(dropnames==Zobs_id)==1){
        themat.Zobs.s = themat; cat("\r","found clique with Zobs!")
        return(themat.Zobs.s)
      }
    }
    new.NEgraph = new.NEgraph[,-drop.ind]
    numleft = dim(new.NEgraph)[2]
    if(length(numleft)==0){ numleft=0 }
    decomp[[ii]] = themat
    cat("\r","found clique",ii,'...')
    ii=ii+1
  }
  return(decomp)
}

#' Generating distance matrix.
#'
#' Generates an example of two-dimensional network of 3 Gaussians with
#' number of units being \code{num_units}.
#' Please refer to the example.
#'
#' @param num_units The number of units in the network.
#'
#' @return A list of \code{D}, a square matrix of all pairwise Euclidean distances between units, and \code{thepoints}, a matrix of dimension \code{num_units} by 2 of coordinates of the generated network units. (binary).
#' @export
out_example_network = function(num_units){

  # generate new network (3 clusters of 2D Gaussians)
  x = c(rnorm(num_units*0.5,sd=0.1)+0.5,rnorm(num_units*0.3,sd=0.075)+0.25,rnorm(num_units*0.2,sd=0.075)+0.3)
  y = c(rnorm(num_units*0.5,sd=0.1)+0.5,rnorm(num_units*0.3,sd=0.075)+0.75,rnorm(num_units*0.2,sd=0.075)+0.3)
  thepoints = cbind(x,y)

  # distance matrix
  eucdist = function(x1,x2){
    sqrt(sum((x1-x2)^2))
  }

  Dmat = array(0,c(nrow(thepoints),nrow(thepoints)))
  for(ii in 1:nrow(thepoints)){
    for(jj in 1:nrow(thepoints)){
      Dmat[ii,jj] = eucdist(thepoints[ii,],thepoints[jj,])
    }
  }
  return(list(D=Dmat,thepoints=thepoints))
}

#' Casts matrix \code{mat} into a sparse matrix.
#'
#' @param mat A matrix to be sparsified
#'
#' @return  A sparse matrix.
#' @export
sparsify <- function(mat){
  mat = Matrix(mat,sparse=TRUE)
  mat
}

#' Calculating average treatment effect
#'
#' Returns difference in means between exposures \code{b} and \code{a}, coded as \code{1} and \code{-1}, respectively.
#'
#' @param Z The treatment vector. Entries of \code{Z} are \code{1} or \code{-1},
#'  indicating that under this treatment, individual is exposed to \code{b} or \code{a}
#' @param Y The outcome vector
#'
#' @return The differenence in means between exposure contrasts \code{b} and \code{a}.
#' @export
ate = function(Z,Y){
  ind1 = which(Z==1)
  ind2 = which(Z==-1)
  mean(Y[ind1])-mean(Y[ind2])
}


#' Generating household structure
#'
#' Define a house structure vector that describes the population in each household,
#' which is useful in the clustered interference.
#' In this running example, houses refer to clusters and kids refer to experimental units.
#'
#' @param N The number of observations (kids)
#' @param K The number of households
#' @param equal A binary parameter that determines the exact house structure. If \code{equal=TRUE},
#'  every household has equal number of kids.
#'  If \code{equal=FALSE}, every household samples the number of kids at random. But in each household,
#'  there would be at least two and at most \code{(2*N/K)-2} kids. Hence, on average there would be \code{N/K} kids per household.
#'
#' @return A vector of length \code{K} that gives the number of kids in each household.
#' @export
out_house_structure = function(N=900,K=N/3,equal=T){
  if(equal){ # exactly N observations
    housestruct = rep(N/K,K)
  }
  if(!equal){ # on average N observations
    housestruct = sample(2:((2*N/K)-2),K,replace=T)
  }
  housestruct
}

#' Get a sample from 1:k
#'
#' @param k the maximum number that we can sample.
#'
#' @return An integer number sampled from 1 to \code{k}.
#' @export
sample_mod = function(k){
  sample(1:k,1)
}

#' Generating exposure status for household
#'
#' Outputs the exposure status for all observations in a household structure. Specifically, given
#' a household structure, it returns which households are to be treated, and within each treated
#' household, which kid is to be treated. Only one kid per treated household is going to be treated.
#'
#' @param housestruct The house structure vector that gives the number of kids in each household
#' @param K1 The number of houses to treat.
#'
#' @return A list specifies treated houses and treated kids in each treated houses.
#' @export
out_treat_household = function(housestruct, K1){
  which_house = sample(1:length(housestruct),size=K1,replace=F) # which households to treat?
  which_kid = apply(as.matrix(housestruct[which_house]),1,sample_mod) # which kids to treat (conditional on the treated households)?
  return(list(which_house=which_house,which_kid=which_kid))
}

#' Gives a range of indices that are exposed to treatments. See \code{out_Z_household}.
#'
#' @param ii the index of interest
#' @param lind the first index in each household
#' @param hind the last index in each household
#'
#' @return A vector of indices exposed to the treatment
#' @export
out_exp_ind = function(ii,lind,hind){
  lb = lind[ii]
  ub = hind[ii]
  lb:ub
}

#' This function returns the treatment status for all observations. Note that treated
#' kid in treated household here also has status \code{1} rather than \code{2}. Please see
#' the example.
#'
#' @param N The number of observations (kids).
#' @param K The number of houses.
#' @param treatment_list A list of treated houses and treated kids.
#' @param housestruct House structure vector that gives number of kids in each household.
#'
#' @return an N x 3 matrix (N = # of observations (kids)), the first column identifies which households are treated with a "1" for the first kid in a treated household, the second and third columns are the ones that matter, they are the treatment and exposure vectors at the observation (kid) level.
#' @export
out_Z_household = function(N,K,treatment_list,housestruct){
  which_house = sort(treatment_list$which_house)
  theord = order(treatment_list$which_house)
  which_kid = treatment_list$which_kid[theord]

  Z = matrix(0,nrow=N,ncol=3)
  colnames(Z) = c('houseind','treat','exp')
  houseind = cumsum(housestruct)-housestruct[1]+1
  Z[houseind,1] = 1 # separates out the different houses
  treatedind = houseind[which_house] + which_kid - 1
  Z[treatedind,2] = 1
  lind = houseind[which_house]
  hind = lind + housestruct[which_house] - 1
  expind = c(apply(as.matrix(1:length(which_house)),1,out_exp_ind,lind=lind,hind=hind))
  Z[expind,3] = 1
  Z
}

#' Returns matrix \code{Z} of dimension (number of units x number of randomizations, i.e. assignments.)
#'
#' @param N The number of observations (kids).
#' @param K The number of houses.
#' @param equal A binary parameter that determines the exact house structure. If equal=TRUE, every household has equal number of kids. If equal=FALSE, every household is constructed by sampling the number of kids at random.
#' @param numrand The number of randomizations
#'
#' @return A matrix of dimension (number of units x number of randomizations, i.e. assignments.) comprised of 0, denoting a unit is pure control under that particular assignment, 1 denoting a unit is a spillover under that particular assignment, and 2, denoting a unit is treated under that particular assignment.
#' @export
out_Zprime = function(N,K,equal=T,numrand){
  housestruct = out_house_structure(N,K,equal)
  Zprime_mat = matrix(0,nrow=N,ncol=numrand)
  for(nn in 1:numrand){
    treatment_list = out_treat_household(housestruct,K1=K/2)
    Zp = out_Z_household(N,K,treatment_list,housestruct)
    Zp_compact = rowSums(Zp[,2:3])
    Zprime_mat[,nn] = Zp_compact # correct the true exposure for those being treated.
  }
  Zprime_mat
}

#' DGP of Basse and Feller (2018)
#'
#' Simulates an outcome vector based on the treatment assignment vector.
#' The specific data generating procedure is adopted from Basse & Feller (2018).
#'
#' @param N The number of observations (kids).
#' @param K The number of houses.
#' @param Zobs Observed treatment assignment vector.
#' @param tau_main The main effect.
#' @param sig_c Standard error on the causal effect. Here \eqn{\sigma_\mu=\sigma_\tau}=\code{sig_c}.
#' @param sig_y Standard error on the outcome vector.
#' @param taus Spillover effect.
#' @param taup Primary effect.
#' @param mu00 The effect on pure control units.
#' @param equal A binary parameter that determines the exact house structure. If \code{equal=TRUE}, every household has equal number of kids. If \code{equal=FALSE}, every household samples the number of kids at random.
#'
#' @return An outcome vector of length \code{N}.
#' @export
out_bassefeller = function(N, K, Zobs, tau_main,
                           sig_c=0.1,
                           sig_y=0.5, taus = 0.7,taup = 1.5,mu00 = 2,equal=T){
  sig_mu <- sig_taup <- sig_taus <- sig_c
  Yi00 = rnorm(K,mu00,sig_mu)
  tauip = rnorm(K,taup,sig_taup)
  tauis = rnorm(K,taus,sig_taus)

  # generate average potential outcomes
  ## Yi11: house i treated and self treated
  ## Yi00: house i not treated and self not treated
  ## Yi10: house i treated and self not treated
  Yi11 =  Yi00 + tauip; Yi10 =  Yi00 + tauis

  # get the experimental structure (OBSERVED)
  housestruct = out_house_structure(N,K,equal)

  # potential outcome vectors
  Yij00 = c()
  Yij10 = c()
  Yij11 = c()
  for(kk in 1:K){
    Yij00 = c(Yij00,rnorm(housestruct[kk],Yi00[kk],sig_y))
    Yij10 = c(Yij10,rnorm(housestruct[kk],Yi10[kk],sig_y))
    Yij11 = c(Yij11,rnorm(housestruct[kk],Yi11[kk],sig_y))
  }

  ## NULL TRUE
  Yij10 = Yij00 + tau_main # enforcing null right here to test validity

  # Yi00 observations
  Yobs = rep(NA, N)
  Yobs = Yij00 * (Zobs==0) + Yij10 * (Zobs==1) + Yij11 * (Zobs==2)
  return(list(Yobs=Yobs,Zobs=Zobs))
}



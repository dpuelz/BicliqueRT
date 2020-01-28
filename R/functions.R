#' One-sided testing
#'
#' Decides to reject or not based on observed test statistic value \code{tobs} and randomization values \code{tvals}.
#'
#' @param tobs The observed value of the test statistic (scalar).
#' @param tvals Vector of randomization values of the test statistic (to compare with \code{tobs}).
#' @param alpha Desired level of the test (between 0 to 1).
#' @param tol Used to check whether \code{tobs} is equal to the 1-\code{alpha} quantile of \code{tvals}.
#' @details
#' The test may randomize to achieve the specified level \code{alpha}
#' when there are very few randomization values.
#' @return Test decision (binary).
#' @seealso Testing Statistical Hypotheses (Ch. 15, Lehman and Romano, 2006)
one_sided_test = function(tobs, tvals, alpha, tol=1e-14) {
  srt = sort(tvals)
  M = length(tvals)
  k = ceiling(M * (1-alpha))
  Tk = srt[k]
  if(abs(tobs - Tk) < tol) {
    # if tobs = Tk
    ax = (M * alpha - sum(tvals > Tk)) / sum(abs(tvals - Tk) < tol)
    return(runif(1) <= ax) ## randomize decision.
  }

  return(tobs > Tk)
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
#'
#' @return Test decision (binary).
#' @seealso Testing Statistical Hypotheses (Ch. 15, Lehman and Romano, 2006)
two_sided_test = function(tobs, tvals, alpha) {
  m1 = one_sided_test(tobs, tvals, alpha=alpha/2)
  m2 = one_sided_test(-tobs, -tvals, alpha=alpha/2) # only one can be 1.
  return(m1 + m2)
}

#' Calculates p-value or test decision
#'
#' Depending on \code{ret_pval} this function returns either a p-value for the test or the binary decision.
#'
#' @param rtest_out A \code{List} with elements \code{tobs}, \code{tvals} (see \link{one_sided_test} for details.)
#' @param ret_pval A \code{Boolean} indicating whether to return a p-value (TRUE) or not.
#' @param alpha Desired test level (from 0 to 1).
#' @return Binary decision if \code{ret_pval} is TRUE, or the p-value otherwise.
#' @details Returns 1 if the test rejects, 0 otherwise.
out_pval = function(rtest_out, ret_pval, alpha) {
  tobs = rtest_out$tobs
  tvals = c(rtest_out$tvals)

  n_all = length(tvals)
  n_higher = sum(tvals > (tobs + 1e-12))
  n_lower = sum(tvals < (tobs - 1e-12))
  n_equal = n_all - n_lower - n_higher

  p1 = (n_equal + n_higher) / n_all  # P(T >= Tobs)
  p2 = (n_equal + n_lower) / n_all  # P(T <= Tobs)

  pval = min(p1, p2)
  if(ret_pval) return(pval)
  return(two_sided_test(tobs, tvals, alpha = alpha))  # this test is less conservative.
}

#' Computes the observed outcome vector, given potential outcomes on the exposures
#'
#' @param z_a A binary vector indicating which units are exposed to \code{a}.
#' @param z_b A binary vector indicating which units are exposed to \code{b}.
#' @param Y_a The potential outcome vector for all units under exposure \code{a}.
#' @param Y_b The potential outcome vector for all units under exposure \code{b}.
#'
#' @return A single observed outcome vector.
out_Yobs = function(z_a,z_b,Y_a,Y_b){
  y = Y_a*z_a + Y_b*z_b
  y
}

#' Returns sparse matrix \code{Z} of dimension (number of units x number of randomizations, i.e. assignments.)
#'
#' @param pi A \code{vector} of length number of units.
#' @param num_randomizations A \code{scalar} denoting the number of assignments desired.
#'
#' @return \code{Z} (sparse binary matrix).  The first column is the observed assignmenht.
out_Z = function(pi,num_randomizations){
  num_units = length(pi)
  Z = Matrix(rbinom(num_randomizations*num_units,1,prob=pi),nrow=num_units,ncol=num_randomizations,sparse=TRUE)
  Z
}

#' One of the main functions for implementing the methodology.  Outputs the null-exposure graph based on binary matrices describing exposure conditions in the null hypothesis.  The null hypothesis is represented as:
#' H_0: Y_i(\code{a}) = Y_i(\code{b}) for all i,
#' and states that potential outcomes are equal for all units exposed to either \code{a} or \code{b}.
#'
#' @param Z_a A binary matrix with dimension (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{a} under assignment j.  Please see example.
#' @param Z_b A binary matrix with (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{b} under assignment j.  Please see example.
#' @param Z A binary matrix of dimension (number of units x number of randomizations, i.e. assignments.) storing the assignment vectors.
#' @param exclude_treated A Boolean denoting whether or not treated units are considered in the hypothesis.  Default is \code{TRUE}.
#'
#' @return \code{NEgraph}, a matrix of dimension (number of units x number of randomizations, i.e. assignments.).  Row i, column j of the matrix either equals -1 (unit i is exposed to \code{a} under assignment j), 1 (unit i is exposed to \code{b} under assignment j), 0 (unit i is neither exposed to \code{a} nor \code{b} under assignment j).
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

#' Identifies the clique to condition on in the randomization test.
#'
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.
#' @param decomp A list containing the clique decomposition of the null-exposure graph.
#'
#' @return The clique to condition on in the randomization test.
out_clique = function(Zobs_id,decomp){
  for(ii in 1:length(decomp)){
    if(sum(Zobs_id==colnames(decomp[[ii]]))){
      return(as.matrix(decomp[ii][[1]]))
    }
  }
}

#' One of the main functions for implementing the methodology. Outputs a clique decomposition of the null-exposure graph.  Specifically, the resulting decomposition paritions the assignment space, while the unit space may overlap.
#'
#' @param NEgraph The null-exposure graph object, see \code{out_NEgraph}.
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.
#' @param minr The minimum number of focal units included in the cliques (algorithm runtime is sensitive to this value).
#' @param minc The minimum number of focal assignment included in the cliques (algorithm runtime is sensitive to this value).
#' @param stop_at_Zobs A Boolean indicating whether the decomposition algorithm should stop when the clique containing the observed assignment is found.  Default value is \code{TRUE}.
#'
#' @return If \code{stop_at_Zobs} is \code{TRUE}, a matrix representing the clique to condition upon. If \code{stop_at_Zobs} is \code{FALSE}, a list containing the clique decomposition of the null-exposure graph.
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
    while(length(themat)==0){
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


#' Generates example two-dimensional network of 3 Gaussians.
#'
#' @param num_units The number of units in the network.
#'
#' @return A list of \code{D}, a square matrix of all pairwise Euclidean distances between units, and \code{thepoints}, a matrix of dimension \code{num_units} by 2 of coordinates of the generated network units. (binary).
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
#' @return  A sparse matrix
sparsify <- function(mat){
  mat = Matrix(mat,sparse=TRUE)
  mat
}

#' Returns difference in means between exposures \code{b} and \code{a}, coded as \code{1} and \code{-1}, respectively.
#'
#' @param Z The treatment vector
#' @param Y The outcome vector
#'
#' @return The differenence in means between exposure contrasts \code{b} and \code{a}
ate = function(Z,Y){
  ind1 = which(Z==1)
  ind2 = which(Z==-1)
  mean(Y[ind1])-mean(Y[ind2])
}

#' The main randomization test function.
#'
#' @param Y The observed outcome vector.
#' @param Z A binary matrix of dimension (number of units x number of randomizations, i.e. assignments.) storing the assignment vectors.
#' @param Z_a A binary matrix with dimension (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{a} under assignment j.  Please see example.
#' @param Z_b A binary matrix with (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{b} under assignment j.  Please see example.
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.
#' @param minr The minimum number of focal units included in the cliques (algorithm runtime is sensitive to this value).
#' @param minc The minimum number of focal assignment included in the cliques (algorithm runtime is sensitive to this value).
#' @param ... Other stuff ...
#'
#' @return A list of items summarizing the randomization test.
#' @examples
#' # generated network - 3 clusters of 2D Gaussians
#' # loads in the 500x500 matrix Dmat (see create_network function ...)
#' # Dmat just encodes all pairwise Euclidean distances between network nodes, and
#' # this is used to define the spillover hypothesis below.
#' set.seed(1)
#' thenetwork = out_example_network(500)
#' D = thenetwork$D
#'
#' # simulation parameters
#' num_randomizations = 5000
#' radius = 0.01
#'
#' # First, construct \code{Z}, \code{Z_a}, \code{Z_b}.
#' # Here, exposure \code{a} is an untreated within \code{radius} of a treated unit, and exposure \code{b} is an untreated unit at least \code{radius} distance away from all treated units.
#' # Experimental design is Bernoulli with prob=0.2.
#' # \code{a_threshold} A scalar denoting the threshold that triggers an exposure to \code{a}.  If exposure \code{a} is simply binary, i.e. whether or not unit j is exposed to \code{a}, then this value should be set to 1.
#' # \code{b_threshold} A scalar denoting the threshold that triggers an exposure to \code{b}.  If exposure \code{b} is simply binary, i.e. whether or not unit j is exposed to \code{b}, then this value should be set to 1.
#' Z = out_Z(pi=rep(0.2,dim(D)[1]),num_randomizations)
#' D_a = sparsify((D<radius)); a_threshold = 1
#' D_b = sparsify((D<radius)); b_threshold = 1
#'
#' Z_a = D_a%*%Z
#' Z_a = sparsify((Z_a>=a_threshold))
#' Z_b = D_b%*%Z
#' Z_b = sparsify((Z_b<b_threshold))
#'
#' simulating an outcome vector
#' Y_a = rnorm(dim(Z)[1])
#' Y_b = Y_a + 0.2
#' Y = out_Yobs(Z_a[,1],Z_b[,1],Y_a,Y_b)
#'
#' run the test
#' CRT = clique_test(Y,Z,Z_a,Z_b,Zobs_id=1,minr=15,minc=15)
clique_test = function(Y,Z,Z_a,Z_b,Zobs_id,minr,minc,...){

  # setting default values if they do not exist
  if(!exists("exclude_treated")){ exclude_treated=TRUE }
  if(!exists("stop_at_Zobs")){ stop_at_Zobs=TRUE }
  if(!exists("ret_pval")){ ret_pval=TRUE }

  # make the null-exposure graph
  cat("construct the null-exposure graph ... \n")
  NEgraph = out_NEgraph(Z_a,Z_b,Z,exclude_treated)

  # decompose the null-exposure graph
  cat("decompose the null-exposure graph ... \n")
  decomp = out_clique_decomposition(NEgraph,Zobs_id,minr,minc,stop_at_Zobs)
  if(stop_at_Zobs){ conditional_clique = decomp }
  if(!stop_at_Zobs){ conditional_clique = out_clique(Zobs_id,decomp) }

  focal_units = as.numeric(rownames(conditional_clique))
  focal_assignments = as.numeric(colnames(conditional_clique))

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

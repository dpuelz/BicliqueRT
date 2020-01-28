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
out_Yobs = function(Z,Y_a,Y_b){
  Za = (Z==-1)
  Zb = (Z==1)
  y = Y_a*Za + Y_b*Zb
  y
}

#' Returns sparse matrix \code{Z} of dimension (number of units x number of randomizations, i.e. assignments.)
#'
#' @param pi A \code{vector} of length number of units.
#' @param num_randomizations A \code{scalar} denoting the number of assignments desired.
#'
#' @return \code{Z} (sparse binary matrix).
out_Z = function(pi,num_randomizations){
  num_units = length(pi)
  Z = Matrix(rbinom(num_randomizations*num_units,1,prob=pi),nrow=samsize,ncol=num_randomizations,sparse=TRUE)
  Z
}

#' One of the main functions for implementing the methodology.  Outputs the null exposure graph based on binary matrices describing exposure conditions in the null hypothesis.  The null hypothesis is represented as:
#' H_0: Y_i(\code{a}) = Y_i(\code{b}) for all i,
#' and states that potential outcomes are equal for all units exposed to either \code{a} or \code{b}.
#'
#' @param D_a A binary square matrix with the number of columns equal to the number of units.  Column i, row j of the matrix corresponds to whether a treated unit j exposes unit i to exposure \code{a}. See examples.
#' @param D_b A binary square matrix with the number of columns equal to the number of units.  Column i, row j of the matrix corresponds to whether a treated unit j exposes unit i to exposure \code{b}. See examples.
#' @param a_threshold A scalar denoting the threshold that triggers an exposure to \code{a}.  If exposure \code{a} is simply binary, i.e. whether or not unit j is exposed to \code{a}, then this value should be set to 1.
#' @param b_threshold A scalar denoting the threshold that triggers an exposure to \code{b}.  If exposure \code{b} is simply binary, i.e. whether or not unit j is exposed to \code{b}, then this value should be set to 1.
#' @param Z A binary matrix of dimension (number of units x number of randomizations, i.e. assignments.) storing the assignment vectors.
#' @param exclude_treated A Boolean denoting whether or not treated units are considered in the hypothesis.  Default is TRUE.
#'
#' @return
out_NEgraph = function(D_a,D_b,a_threshold,b_threshold,Z,exclude_treated=TRUE){
  # first, compute Z_a matrix

  # second, compute Z_b matrix

  # if exclude_treated==TRUE, then component multiply D_a, D_b by !Z.

  # construct NEgraph from the exposures matrices

  return(NEgraph)
}

out_clique = function(z.id,decom){
  for(ii in 1:length(decom)){
    if(sum(z.id==colnames(decom[[ii]]))){
      return(ii)
    }
  }
}

out_clique_decomposition = function(NEgraph,Zobs_id,minr,minc,stop_at_Zobs=FALSE){

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
        themat.Zobs.s = themat; cat("found clique with Zobs!\n")
        return(themat.Zobs.s)
      }
    }
    new.NEgraph = new.NEgraph[,-drop.ind]
    numleft = dim(new.NEgraph)[2]
    if(length(numleft)==0){ numleft=0 }
    decomp[[ii]] = themat
    cat("found clique",ii,'...\n')
    ii=ii+1
  }
  return(decomp)
}

clique_test = function(Yobs,BImat,Zobs_id){
  focalunits = as.numeric(rownames(BImat))
  cliqassign = as.numeric(colnames(BImat))

  # restrict Y to clique
  Y.clique = Yobs[focalunits]

  # test statistics
  Zobs_cliqid = which(Zobs_id==cliqassign)
  Tobs = ate(BImat[, Zobs_cliqid], Y.clique)

  Trand = apply(as.matrix(BImat[, -Zobs_cliqid]), 2, function(z) ate(z, Y.clique))

  reject = (one_sided_test(tobs=Tobs,tvals=Trand,alpha=0.05,tol=1e-14))
  reject
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
sparsify <- function(mat){
  mat = Matrix(mat,sparse=TRUE)
  mat
}

out_DZ = function(r,D,Z){
  Dr = binarize_r(r,D)
  Z = sparsify(Z)
  DZ = Dr%*%Z
  DZ
}

out_binaryDZ = function(DZ,direction){
  if(direction=='lt'){
    DZbin = (DZ>0)
  }
  if(direction=='gt'){
    DZbin = (DZ==0)
  }
  DZbin
}

binarize_r <- function(r,mat){
  # r is a distance metric, such as euclidean distance between street segments
  # every segment pair > r distance apartment will be zero.
  matnew = Matrix(0,dim(mat)[1],dim(mat)[2],sparse=TRUE)
  rownames(matnew) = rownames(mat)
  colnames(matnew) = colnames(mat)

  matnew[mat<=r] <- 1
  diag(matnew) <- -1e10 # this is key here .. focus hypothesis on untreated units only
  matnew
}

#' Functions that returns difference in means between exposures \code{b} and \code{a}, coded as \code{1} and \code{-1}, respectively.
ate = function(Z,Y){
  ind1 = which(Z==1)
  ind2 = which(Z==-1)
  mean(Y[ind1])-mean(Y[ind2])
}

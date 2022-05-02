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
#' @param tau_equal If \code{TRUE}, then regardless of structual parameters, we set Y_i(1) = Y_i(0) + \code{tau_main}.
#' @param tau_main The main effect such that we force Y_i(1) = Y_i(0) + \code{tau_main}.
#' @param sig_c Standard error on the causal effect. Here \eqn{\sigma_\mu=\sigma_\tau}=\code{sig_c}.
#' @param sig_y Standard error on the outcome vector.
#' @param taus Spillover effect.
#' @param taup Primary effect.
#' @param mu00 The effect on pure control units.
#' @param equal A binary parameter that determines the exact house structure. If \code{equal=TRUE}, every household has equal number of kids. If \code{equal=FALSE}, every household samples the number of kids at random.
#'
#' @return An outcome vector of length \code{N}.
#' @export
out_bassefeller = function(N, K, Zobs, tau_main, tau_equal = F,
                           sig_c = 0.1, sig_y = 0.5, taus = 0.7, taup = 1.5, mu00 = 2,
                           equal = T){
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
  if (tau_equal){
    Yij10 = Yij00 + tau_main
  }

  # Yi00 observations
  Yobs = rep(NA, N)
  Yobs = Yij00 * (Zobs==0) + Yij10 * (Zobs==1) + Yij11 * (Zobs==2)
  return(list(Yobs=Yobs,Zobs=Zobs))
}


#' Diagnostic for Randomized Experiment
#'
#' For different number of units selected, it returns an \code{fmax} by 2 matrix,
#' where \code{fmax} is the maximum number of units to be selected (see below).
#' The first column is the approximate number of focal assignments of the largest possible biclique
#' of the null exposure graph that takes the selected units as focal units if we draw
#' \code{N} randomizations. The second column is the first column with an additional requirement
#' that each of the biclique's focal assignments does not contain only one type of exposure.
#' It can serve as a diagnostic for choosing \code{minr}, \code{minc} or \code{minass} in the
#' function \code{clique_test}.
#'
#' @param struc Structure of the cluster. It should have two columns. The first column
#' represents group ID, and the second column is the individual ID. Please see the example.
#' @param p_group The probability of a group being selected as treated group.
#' @param p_indi_t The probability of an individual in the treated group being treated.
#' @param p_indi_nt The probability of an individual in the non-treated group being treated.
#' @param N The number of randomizations performed.
#' @param fmax Maximum number of units to be selected to do the diagnosis. Default is 25.
#' @return an \code{fmax} by 2 matrix.
#' @export
clique_diagnostic = function(struc, p_group, p_indi_t, p_indi_nt, N, fmax = 25, NR = 500, Nx = 500){
  K = length(unique(struc[,1]))
  p_a = (1 - p_indi_t); p_b = (1 - p_indi_nt)
  prob_group = tapply(struc[,2], factor(struc[,1]), length) / dim(struc)[1]
  largest_clique = matrix(nrow = fmax, ncol = 2)

  for (f in 1:fmax){
    cat('\r','doing diagnostic for the number of focal units being ',f,'...')
    prob_f = 0; prob_f_val = 0
    x_rand = rmultinom(n = Nx, size = f, prob_group)

    for (idx in 1:Nx){
      x_rand_idx = x_rand[,idx]
      x_rand_idx_selected = which(x_rand_idx==1)
      R_idx = replicate(NR, rbinom(n = K, size = 1, 0.5), simplify = FALSE)
      FR_idx = 0; FR_idx_val = 0
      for (idr in 1:NR){
        FR_idr = (p_a^x_rand_idx - p_b^x_rand_idx) * R_idx[[idr]] + p_b^x_rand_idx
        R_selected = R_idx[[idr]][x_rand_idx_selected]

        FR_idx = FR_idx + prod(FR_idr)
        FR_idx_val = FR_idx_val + prod(FR_idr) * ((sum(R_selected)!=f) & (sum(R_selected)!=0))
      }
      FR_idx = FR_idx / NR; FR_idx_val = FR_idx_val / NR
      prob_f = prob_f + FR_idx; prob_f_val = prob_f_val + FR_idx_val
    }
    prob_f = prob_f / Nx; prob_f_val = prob_f_val / Nx
    largest_clique[f, 1] = prob_f * N
    largest_clique[f, 2] = prob_f_val * N
  }

  return(largest_clique)
}

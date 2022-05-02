#' The main randomization test function for contrast hypothesis
#'
#' @param Y The observed outcome vector.
#' @param Z A binary matrix of dimension (number of units x number of randomizations, i.e. assignments.) storing the assignment vectors. Please see example.
#' @param Z_a A binary matrix with dimension (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{a} under assignment j. Please see example.
#' @param Z_b A binary matrix with (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{b} under assignment j. Please see example.
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.
#' @param Xadj The covariates that might affect Y. If not \code{NULL}, will replace \code{Y} by the residuals from the linear regression of \code{Y} on \code{Xadj}. Note that users would need to add an intercept to \code{Xadj} manually if they want.
#' To adjust \code{Xadj}, pass in \code{Y_adj=TRUE} and a non-empty \code{NULL} that has the same row numbers as \code{Y}.
#' @param alpha The significance level. By default it's \eqn{0.05}.
#' @param tau The \eqn{\tau} in the null \eqn{Y_i(b) = Y_i(a) + \tau} for all \eqn{i}. By default is \eqn{0}.
#' @param decom The algorithm used to calculate the biclique decomposition. Currently supported algorithms
#' are "bimax" and "greedy".
#' @param ret_ci Whether calculates the \eqn{1-\alpha} confidence interval or not. Default is \code{FALSE}.
#' @param ... Other stuff ...
#'
#' @return A list of items summarizing the randomization test. If for some focal assignments
#' in the biclique that contains \code{Zobs}, exposures for each unit are the same, it will
#' contain an error message, and the test decision will be \code{NA}.
#' @export
clique_test = function(Y, Z, Z_a, Z_b, Zobs_id, Xadj=NULL, alpha=0.05, tau=0,
                       decom='bimax', ret_ci=FALSE, ci_dec=2, ci_method='grid', ...){
# clique_test = function(Y, Z, Z_a, Z_b, Zobs_id, Xadj=NULL, alpha=0.05, tau=0,
#                        decom='bimax', ret_ci=FALSE, ci_dec=2, ci_method='grid', ...){
  addparam = list(...) # catch variable parameters

  # setting default values if they do not exist
  if(is.null(addparam$exclude_treated)){ exclude_treated=TRUE } else { exclude_treated=addparam$exclude_treated }
  if(is.null(addparam$stop_at_Zobs)){ stop_at_Zobs=FALSE } else { stop_at_Zobs= addparam$stop_at_Zobs}
  if(is.null(addparam$ret_pval)){ ret_pval=TRUE } else { ret_pval=addparam$ret_pval }
  if(is.null(addparam$adj_Y)){ adj_Y=FALSE } else { adj_Y=(addparam$adj_Y)&(!is.null(Xadj)) } # is TRUE iff we specify it and Xadj is not null

  # make the null-exposure graph
  cat("construct the null-exposure graph ... \n")
  NEgraph = out_NEgraph(Z_a,Z_b,Z,exclude_treated)

  # decompose the null-exposure graph
  cat("decompose the null-exposure graph ... \n")
  if (decom == 'bimax'){
    decomp = out_clique_decomposition_bimax(NEgraph, Zobs_id,
                                            minr=addparam$minr, minc=addparam$minc, stop_at_Zobs)
    if(stop_at_Zobs){
      conditional_clique = decomp
    } else {conditional_clique = out_clique(Zobs_id,decomp)}
  }

  if (decom == 'greedy'){
    decomp = out_clique_decomposition_greedy(NEgraph, Zobs_id, num_ass=addparam$minass, stop_at_Zobs)
    if(stop_at_Zobs){
      conditional_clique = decomp
    } else {conditional_clique = out_clique(Zobs_id,decomp)}
  }

  focal_units = as.numeric(rownames(conditional_clique)) # a list of row names
  focal_assignments = as.numeric(colnames(conditional_clique)) # a list of column names


  # adjust for covariates, if any
  # currently can only be adjusted by linear regression
  if (adj_Y){
    Y = c(Y - Xadj %*% solve(t(Xadj) %*% Xadj) %*% t(Xadj) %*% Y)
  }

  # check whether the conditional_clique has non-unique exposure for each focal assignment, which is essential for the test.
  check_clique_unique <- apply(conditional_clique, 2, function(x) length(unique(x))==1)
  if (sum(check_clique_unique)!=0){
    cat("The biclique containing Zobs contains only one exposure for some focal assignment, and therefore the result is a powerless test \n")
    check_clique_unique <- "Powerless test, because the biclique containing Zobs contains only one exposure for some focal assignment"
  } else {
    check_clique_unique <- c()
  }

  if (!ret_ci){
    # test the null hypothesis

    # adjust Y for null hypothesis
    Y = out_Yobs(Z_a, Z_b, Y, Y+tau)

    # restrict Y to clique
    Y.clique = Y[focal_units]

    # run the test
    cat("\n")
    cat("run the clique-based randomization test ... \n")
    Zobs_cliqid = which(Zobs_id==focal_assignments)
    tobs = ate(conditional_clique[, Zobs_cliqid], Y.clique)
    tvals = apply(as.matrix(conditional_clique[, -Zobs_cliqid]), 2, function(z) ate(z, Y.clique))
    rtest_out = list(tobs=tobs, tvals=tvals)
    decision = out_pval(rtest_out, ret_pval, alpha)

    if (sum(abs(tvals-tobs))/length(focal_assignments) < 1e-14){
      cat("The test statistics are identical under each focal assignment. p-value is not available")
      check_equal_statistics <- "The test statistics are identical under each focal assignment. p-value is not available"
      decision = NA
    } else {
      check_equal_statistics <- c()
    }

    # organize into return list
    retlist = list(decision=decision, ret_pval=ret_pval, tobs=tobs,
                   tvals=tvals, focal_units=focal_units, focal_assignments=focal_assignments,
                   NEgraph=NEgraph, warnings=c(check_clique_unique, check_equal_statistics))
    retlist
  } else {
    # return a confidence interval by inverting the test

    Zobs_cliqid = which(Zobs_id==focal_assignments)
    if (is.null(addparam$ci_ub)) {ci_ub = as.numeric(quantile(Y[Z_b[,Zobs_id]], 0.95)-quantile(Y[Z_a[,Zobs_id]], 0.05))} else {ci_ub = addparam$ci_ub}
    if (is.null(addparam$ci_lb)) {ci_lb = as.numeric(quantile(Y[Z_b[,Zobs_id]], 0.05)-quantile(Y[Z_a[,Zobs_id]], 0.95))} else {ci_lb = addparam$ci_lb}
    ci_ub = round(ci_ub, ci_dec) + 10^(-ci_dec+1) # ci_dec is an integer specifying decimals for CI. By default is 2.
    ci_lb = round(ci_lb, ci_dec) - 10^(-ci_dec+1)

    if (ci_method == "grid"){
      ci_out = CI_grid(ci_lb, ci_ub, ci_dec)
    } else if (ci_method == "bisection"){
      ci_out = CI_bisection(ci_lb, ci_ub, ci_dec)
    }

    ci_lb_out = ci_out[1]
    ci_ub_out = ci_out[2]
    cat(sprintf("The %s%% confidence interval is [%s,%s] \n", round((1-alpha)*100, 4), ci_lb_out, ci_ub_out))

    retlist = list(ci=list(lb=ci_lb_out, ub=ci_ub_out),
                   focal_units=focal_units, focal_assignments=focal_assignments,
                   NEgraph=NEgraph, warnings=check_clique_unique)
    retlist
  }

}



#' The main randomization test function for the full exclusion null hypothesis.
#'
#' @param Y The observed outcome vector.
#' @param Zrealized A vector of length being number of units that gives the realization of treatment assignment.
#' @param design_fn A function whose return is a one time realization of the distribution of the treatment.
#' @param exposure_fn A function that specifies the exposure of each individual under different randomization.
#' @param decom The algorithm used to calculate the biclique decomposition. Currently supported algorithms are "bimax" and "greedy".
#' @param num_randomizations Number of randomizations to perform.
#' @param ... Other stuff, such as parameters for clique decomposition algorithms (\code{minass} for
#' greedy decomposition, \code{minr} and \code{minc} for bimax decomposition).
#'
#' @return A list of items summarizing the randomization test.
clique_test_ex = function(Y, Zrealized, design_fn, exposure_fn, decom="greedy", alpha=0.05,
                          num_randomizations=2000, ...){

  # design_fn replaces "Z". When generating realizations, set the ture realization be 1: Zobs_id=1
  # num_rand is given currently, in the future perhaps automatically select.
  # exposure_fn replaces "expos"
  # we now give directly Zrealized: #unit x 1 vector, so no need Zobs_id


  addparam = list(...) # catch variable parameters

  # setting default values if they do not exist
  if(is.null(addparam$exclude_treated)){ exclude_treated=TRUE } else { exclude_treated=addparam$exclude_treated }
  if(is.null(addparam$stop_at_Zobs)){ stop_at_Zobs=FALSE } else { stop_at_Zobs= addparam$stop_at_Zobs}
  if(is.null(addparam$ret_pval)){ ret_pval=TRUE } else { ret_pval=addparam$ret_pval }
  if(is.null(addparam$adj_Y)){ adj_Y=FALSE } else { adj_Y=(addparam$adj_Y)&(!is.null(Xadj)) } # is TRUE iff we specify it and Xadj is not null

  # generate randomizations using design_fn & num_randomizations --> Z
  Z = matrix(0, nrow=length(Zrealized), ncol=(num_randomizations+1))
  Z[, 1] = Zrealized
  for (id_rand in 1:num_randomizations){
    Z[, id_rand+1] = design_fn()
  }
  Zobs_id = 1

  # generate exposure using exposure_fn --> expos
  num_units = length(Y)
  dim_exposure = length(exposure_fn(Zrealized, 1)) # get how long the exposure is

  if (dim_exposure > 1){
    expos = array(0, c(num_units, num_randomizations+1, dim_exposure))
  } else if (dim_exposure == 1){
    expos = array(0, c(num_units, num_randomizations+1, 2)) # to make sure expos is always 3-D
  } else {
    cat("exposure_fn should output a vector of length no less than 1!")
  }
  for (id_rand in 1:(num_randomizations+1)){
    for (id_unit in 1:num_units){
      expos[id_unit, id_rand, ] = exposure_fn(Z[,id_rand], id_unit)
    }
  }

  # decompose the null-exposure graph
  cat("decompose the null-exposure graph ... \n")
  Z0 = 1:dim(Z)[2]
  conditional_clique_idx = list(focal_unit=c(), focal_ass=c())
  conditional_clique = Matrix(0, nrow=dim(Z)[1], ncol=dim(Z)[2], sparse=TRUE)

  # while length(Z0)>0, do:
  # 1. select one random from 1:length(Z0) as Zsub
  # 2. generate multiNEgraph = out_NEgraph_ex(Zsub, Z0, expos)
  # 3. decompose multiNEgraph to multi_clique, get one biclique is enough
  # 4. conditional_clique =union of multi_clique, delete multi_clique's col from Z0
  if (decom == 'bimax'){
    while (length(Z0)>0){
      cat("\r","Z0 now has length", length(Z0), '...')

      Zsub = sample(1:length(Z0), size = 1) # Zsub here is an index of vector Z0
      multiNEgraph = out_NEgraph_ex(Zsub, Z0, expos)

      iremove = which(rowSums(multiNEgraph!=0)==0)  # removes isolated units.
      if(length(iremove)!=0){ multiNEgraph = multiNEgraph[-iremove,] }

      numleft = ncol(multiNEgraph)
      minr.new = min(addparam$minr, numleft)
      minc.new = min(addparam$minc, numleft)
      bitest = biclust(multiNEgraph, method=BCBimax(), minr.new, minc.new, number=1)
      bicliqMat = bicluster(multiNEgraph, bitest)
      themat = bicliqMat$Bicluster1

      # if current minr & minc gives no biclique decom, try a smaller one
      while((length(themat)==0) & (numleft > 1)){
        numleft = numleft - 1
        minr.new = min(addparam$minr, numleft)
        minc.new = min(addparam$minc, numleft)
        bitest = biclust(multiNEgraph, method=BCBimax(), minr.new, minc.new, number=1)
        bicliqMat = bicluster(multiNEgraph, bitest)
        themat = bicliqMat$Bicluster1
      }

      if (is.matrix(themat)){
        focal_unit = as.integer(rownames(themat))
        focal_ass = as.integer(colnames(themat)); focal_ass_match = Z0[focal_ass]
      } else { # the biclique we get is only one single column where Zsub lies (b/c this col is all 1)
        next # perhaps just skip this loop and redraw a new Zsub
      }
      conditional_clique[focal_unit, focal_ass_match] = 1
      conditional_clique_idx$focal_unit = union(conditional_clique_idx$focal_unit, focal_unit)
      conditional_clique_idx$focal_ass = union(conditional_clique_idx$focal_ass, focal_ass_match)
      Z0 = Z0[-focal_ass]

      # stop at Zobs: can stop when one of the biclique we found contains Zobs
      stop_at_Zobs = T
      if (stop_at_Zobs){
        if (sum(focal_ass_match==Zobs_id)>0){
          break
        }
      }
    }
  }

  if (decom == 'greedy'){
    while (length(Z0)>0){
      cat("\r","Z0 now has length", length(Z0), '...')

      failed_Zsub = c() # record Zsub that cannot give a greedy clique, to speed up the decomposition a bit
      while (length(failed_Zsub) < length(Z0)){
        Zsub = sample(setdiff(1:length(Z0), failed_Zsub), size = 1) # Zsub here is an index of vector Z0
        multiNEgraph = out_NEgraph_ex(Zsub, Z0, expos)
        num_ass = addparam$minass
        break_signal = FALSE

        iremove = which(rowSums(multiNEgraph!=0)==0)  # removes isolated units.
        if(length(iremove)!=0){ multiNEgraph = multiNEgraph[-iremove,] }

        if (dim(multiNEgraph)[2]<=num_ass){ # when remaining cols not enough to do the greedy algo.
          units_leftover = which(rowSums(multiNEgraph^2)==dim(multiNEgraph^2)[2])
          themat = multiNEgraph[units_leftover,]
          break_signal = TRUE
        } else {
          test = out_greedy_decom(multiNEgraph, num_ass)
          themat = test$clique
        }
        if (is.matrix(themat)) {
          if (dim(themat)[1]>0){
            focal_unit = as.integer(rownames(themat))
            focal_ass = as.integer(colnames(themat)); focal_ass_match = Z0[focal_ass]
            break # break the (length(failed_Zsub) < length(Z0)) loop
          } else { # the decomposed clique is a 0xn matrix
            failed_Zsub = c(failed_Zsub, Zsub)
            next
          }
        } else { # the decomposed clique is not a matrix, skip to next loop
          failed_Zsub = c(failed_Zsub, Zsub)
          next
        }
      }

      conditional_clique[focal_unit, focal_ass_match] = 1
      conditional_clique_idx$focal_unit = union(conditional_clique_idx$focal_unit, focal_unit)
      conditional_clique_idx$focal_ass = union(conditional_clique_idx$focal_ass, focal_ass_match)
      Z0 = Z0[-focal_ass]

      stop_at_Zobs = TRUE
      if (stop_at_Zobs){
        if (sum(focal_ass_match==Zobs_id)>0){
          break
        }
      }
      if (break_signal) {break}
    }
  }

  # test
  cat("\n finding test statistics ... \n")
  conditional_clique = as.matrix(conditional_clique)
  teststat = apply(conditional_clique, 2, function(z) {mean(Y[as.logical(z)])})
  tobs = teststat[Zobs_id]; tvals = teststat[setdiff(conditional_clique_idx$focal_ass, Zobs_id)]
  pval = out_pval(list(tobs=tobs, tvals=tvals), T, alpha) # here we use the previous out_pval function

  retlist = list(pval=pval, tobs=tobs, tvals=tvals, conditional_clique=conditional_clique)
  return(retlist)

}




#' The generalized main randomization test function.
#'
#' @param Yrealized The observed outcome vector.
#' @param Zrealized A vector that gives the realization of treatment assignment. Its length should match \code{Yrealized}
#' and is the number of units in the experiment.
#' @param hypothesis A list that contains three functions specifyting the experiment design and null hypothesis.
#' See details for further illustration.
#' @param controls A list that contains settings for biclique decomposition and covariates adjustment.
#' By default it is \code{list(method="greedy", minass=10, num_randomizations=2000)}.
#' That is, by default it uses greedy decomposition algorithm with \code{minass=10}, and the number of
#' randomizations to perform is 2000. See details for further illustration.
#' @param teststat The test statistic used. It should contain at least (with order)
#' \code{Y, Z, focal_unit_indicator} as inputs where \code{Y} is the outcome vector, \code{Z} is the treatment vector
#' and \code{focal_unit_indicator} is a 0-1 vector indicating whether a unit is focal (=1) or not (=0). All three inputs
#' should have length equal to number of units and have the same ordering.
#' @param alpha The significance level. By default it's \code{0.05}.
#'
#' @details \code{hypothesis} contains three functions:
#' \itemize{
#'  \item{\code{design_fn}} {A function that returns a realization of treatment for the whole sample. For example,
#' if each unit has equal probability \eqn{0.2} to receive the treatment independently, we can write
#' \code{design_fn = function() { rbinom(num_units, 1, prob=0.2) }}.}
#'  \item{\code{exposure_i}} {A function that returns exposure \eqn{f_i(z)} of unit \eqn{i} under treatment \eqn{z} where
#' \eqn{z} is the treatment for the whole sample. The inputs of the function are an index \code{i} and
#' a vector \code{z}. For example, if the exposure of \eqn{i} under \eqn{z} is the treatment it receives, then
#' we can write \code{exposure_i = function(z, i) { z[i] }}. See more examples in the README file.}
#'  \item{\code{null_equiv}} {A function that takes two inputs from \code{exposure_i} and determines
#' whether \eqn{f_i(z_1)} is equivalent to \eqn{f_i(z_2)} under the null hypothesis. For example, if the
#' null is "extent of interference" type of null, we can write
#' \code{null_equiv = function(e1, e2) {identical(e1, e2)}}.}
#' }
#' \code{controls} contains several components:
#' \itemize{
#'  \item{\code{method}} {Specifies the decomposition method. Should be either \code{"bimax"} or \code{"greedy"}.}
#'  \item{\code{minr},\code{minc} or \code{minass}} {If \code{"bimax"} is used, \code{minr} and \code{minc}
#'  should be supplied that specify the minimum number of units and assignments
#'  in the bicliques found by the algorithm. If \code{"greedy"} is used, \code{minass} should be supplied.}
#'  \item{\code{num_randomizations}} {Number of randomizations to perform. If it is not specified, will be set
#'  to be 2000 by default.}
#'  \item{(optional) \code{Xadj}} {The covariates that might affect Y. If it is specified in \code{controls},
#'   will replace \code{Yrealized} by the residuals from the linear regression of \code{Yrealized} on \code{Xadj}
#'   (number of rows in \code{Yrealized} and in \code{Xadj} should be the same).
#'   Note that users would need to add an intercept to \code{Xadj} manually if they want.}
#' }
#'
#' @return A list of items summarizing the randomization test. It contains the p-value \code{p.value},
#' test statistic \code{statistic}, the randomization distribution of the test statistic \code{statistic.dist},
#' and a list \code{MNE} of clique decomposition. Each element of \code{MNE} records one biclique decomposed
#' from a multi-null exposure graph and contains its focal units and focal assignments.
#'
clique_test = function(Yrealized, Zrealized, hypothesis, teststat, alpha=0.05,
                       controls=list(method="greedy", minass=10, num_randomizations=2000)){

  # catch functions and controls from hypothesis and controls
  design_fn = hypothesis$design_fn
  exposure_i = hypothesis$exposure_i
  null_equiv = hypothesis$null_equiv
  num_randomizations = controls$num_randomizations
  if (is.null(num_randomizations)) {num_randomizations = 2000} # set to be 2000 if not supplied

  # check length of Yrealized and Zrealized
  stopifnot("length of Yrealized and Zrealized should be the same" = (length(Yrealized) == length(Zrealized)))

  # check components in hypothesis
  stopifnot("hypothesis should contain all of design_fn, exposure_i and null_equiv"
            = (!(is.null(design_fn) | is.null(exposure_i) | is.null(null_equiv))) )

  # check decomposition in controls
  decom = controls$method
  if (is.null(decom)) {
    stop("controls should contain the decomposition method", call. = F)
  } else if (decom=="bimax") {
    minr = controls$minr
    minc = controls$minc
    if (!(is.numeric(minr)&is.numeric(minc))) {
      stop("if 'bimax' is used in controls, should supply both 'minr' and 'minc' as integer", call. = F)
    }
  } else if (decom=="greedy"){
    minass = controls$minass
    if (!is.numeric(minass)) {
      stop("if 'greedy' is used in controls, should supply 'minass' as integer", call. = F)
    }
  } else {
    stop("the decomposition method in controls should be either 'bimax' or 'greedy'", call. = F)
  }

  # generate randomizations using design_fn & num_randomizations --> Z_m
  num_units = length(Yrealized)
  Z_m = matrix(0, nrow=num_units, ncol=(num_randomizations+1))
  Z_m[, 1] = Zrealized
  for (id_rand in 1:num_randomizations){
    Z_m[, id_rand+1] = design_fn()
  }
  Zobs_id = 1

  # adjust for covariates, if any. Currently can only be adjusted by linear regression
  Xadj = controls$Xadj
  if (!is.null(Xadj)){
    Yrealized = c(Yrealized - Xadj %*% solve(t(Xadj) %*% Xadj) %*% t(Xadj) %*% Yrealized)
  }

  # generate exposure for each unit under different treatment assignment
  dim_exposure = length(exposure_i(Zrealized, 1)) # get how long the exposure is
  if (dim_exposure > 1){
    expos = array(0, c(num_units, num_randomizations+1, dim_exposure))
  } else if (dim_exposure == 1){
    expos = array(0, c(num_units, num_randomizations+1, 2)) # to make sure expos is always 3-D
  } else {
    stop("exposure_i should output a vector of length no less than 1")
  }
  for (id_rand in 1:(num_randomizations+1)){
    for (id_unit in 1:num_units){
      expos[id_unit, id_rand, ] = exposure_i(Z_m[,id_rand], id_unit)
    }
  }

  # decompose the null-exposure graph
  cat("decompose the null-exposure graph ... \n")
  Z0 = 1:dim(Z_m)[2]
  MNE = list()

  # while length(Z0)>0, do:
  # 1. select one random from 1:length(Z0) as Zsub
  # 2. generate multiNEgraph = out_NEgraph_multi(Zsub, Z0, expos)
  # 3. decompose multiNEgraph to multi_clique, get one biclique is enough
  # 4. conditional_clique =union of multi_clique, delete multi_clique's col from Z0
  if (decom == 'bimax'){
    method = "Clique test with Bimax decomposition."
    while (length(Z0)>0){
      cat("\r","Z0 now has length", length(Z0), '...')

      Zsub = sample(1:length(Z0), size = 1) # Zsub here is an index of vector Z0
      # multiNEgraph = out_NEgraph_multi(Zsub, Z0, Z_m, num_units, exposure_fn)
      multiNEgraph = out_NEgraph_multi_separate(Zsub, Z0, expos, null_equiv, dim_exposure)

      iremove = which(rowSums(multiNEgraph!=0)==0)  # removes isolated units.
      if(length(iremove)!=0){ multiNEgraph = multiNEgraph[-iremove,] }

      numleft = ncol(multiNEgraph)
      minr.new = min(minr, numleft)
      minc.new = min(minc, numleft)
      bitest = biclust(multiNEgraph, method=BCBimax(), minr.new, minc.new, number=1)
      bicliqMat = bicluster(multiNEgraph, bitest)
      themat = bicliqMat$Bicluster1

      # if current minr & minc gives no biclique decom, try a smaller one
      while((length(themat)==0) & (numleft > 1)){
        numleft = numleft - 1
        minr.new = min(minr, numleft)
        minc.new = min(minc, numleft)
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

      Z_m_assignments = Z_m[,focal_ass_match]
      rownames(Z_m_assignments) = 1:num_units
      colnames(Z_m_assignments) = focal_ass_match
      MNE = append(MNE, list(list(units = focal_unit, assignments = Z_m_assignments)))
      Z0 = Z0[-focal_ass]

      # stop when one of the biclique we found contains Zobs
      if (sum(focal_ass_match==Zobs_id)>0){
        break
      }
    }
  }

  if (decom == 'greedy'){
    method = "Clique test with greedy decomposition."
    while (length(Z0)>0){
      cat("\r","Z0 now has length", length(Z0), '...')

      failed_Zsub = c() # record Zsub that cannot give a greedy clique, to speed up the decomposition a bit
      while (length(failed_Zsub) < length(Z0)){
        Zsub = sample(setdiff(1:length(Z0), failed_Zsub), size = 1) # Zsub here is an index of vector Z0
        # multiNEgraph = out_NEgraph_multi(Zsub, Z0, Z_m, num_units, exposure_fn)
        multiNEgraph = out_NEgraph_multi_separate(Zsub, Z0, expos, null_equiv, dim_exposure)
        num_ass = minass
        break_signal = FALSE

        ### should not remove isolated units here, otherwise rownames of multiNEgraph is distorted,
        ### then get_clique function will go wrong when matching rownames.
        # iremove = which(rowSums(multiNEgraph!=0)==0)  # removes isolated units.
        # if(length(iremove)!=0){ multiNEgraph = multiNEgraph[-iremove,] }

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

      Z_m_assignments = Z_m[,focal_ass_match]
      rownames(Z_m_assignments) = 1:num_units
      colnames(Z_m_assignments) = focal_ass_match
      MNE = append(MNE, list(list(units = focal_unit, assignments = Z_m_assignments)))
      Z0 = Z0[-focal_ass]

      # stop when one of the biclique we found contains Zobs
      if (sum(focal_ass_match==Zobs_id)>0){
        break
      }
      if (break_signal) {break}
    }
  }

  # test, focal_unit and focal_ass_match
  cat("\n finding test statistics ... \n")
  focal_unit_indicator = 1:num_units %in% focal_unit # indicator of whether a unit is focal
  test_stats = rep(0, length(focal_ass_match))
  for (zid in 1:length(focal_ass_match)){
    Z_focal = Z_m[,focal_ass_match[zid]]
    test_stats[zid] = teststat(Yrealized, Z_focal, focal_unit_indicator)
  }
  tobs = test_stats[which(focal_ass_match==Zobs_id)]
  tvals = test_stats[which(focal_ass_match!=Zobs_id)]
  pval = out_pval(list(tobs=tobs, tvals=tvals), T, alpha) # here we use the previous out_pval function

  if (sum(abs(tvals-tobs))/length(test_stats) < 1e-14){
    warning("Randomization distribution is degenerate. p.value is not available. Consider changing the test statistic.
            See ... for further discussions.")
  }

  # return
  retlist = list(p.value=pval, statistic=tobs, statistic.dist=tvals, method=method, MNE=MNE)
  return(retlist)

}






#' The main randomization test function for contrast hypothesis.
#'
#' @param Y The observed outcome vector.
#' @param Z A binary matrix of dimension (number of units x number of randomizations, i.e. assignments.) storing the assignment vectors. Please see example.
#' @param Z_a A binary matrix with dimension (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{a} under assignment j. Please see example.
#' @param Z_b A binary matrix with (number of units x number of randomizations, i.e. assignments.)  Row i, column j of the matrix corresponds to whether a unit i is exposed to \code{b} under assignment j. Please see example.
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.
#' @param Xadj The covariates that might affect Y. If not \code{NULL}, will replace \code{Y} by the residuals from the linear regression of \code{Y} on \code{Xadj}. Note that users would need to add an intercept to \code{Xadj} manually if they want.
#' To adjust \code{Xadj}, pass in \code{adj_Y=TRUE} and a non-empty \code{NULL} that has the same row numbers as \code{Y}.
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
clique_test_contrast = function(Y, Z, Z_a, Z_b, Zobs_id, Xadj=NULL, alpha=0.05, tau=0,
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
  NEgraph = out_NEgraph_contrast(Z_a,Z_b,Z,exclude_treated)

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



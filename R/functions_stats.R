#' One-sided testing
#'
#' Compute a one-sided p-value given observed test statistic \code{tobs}
#' and randomized \code{tvals}.
#'
#' @param tobs The observed value of the test statistic (scalar).
#' @param tvals A vector of randomization values of the test statistic (to compare with \code{tobs}).
#'
#' @details The test may randomize the p-value depending on how many \code{tvals}
#'  are the same as \code{tobs} (up to a tolerance level) to get better performance.
#'
#' @return p-value for one-sided testing.
#' @seealso Testing Statistical Hypotheses (Ch. 15, Lehman and Romano, 2006)
#' @export
one_sided_test = function(tobs, tvals, tol=1e-14) {
  M = sum(abs(tvals-tobs)<tol)
  p_val = (sum(tvals>tobs) + runif(1)*(1+M)) / (length(tvals) + 1)
  p_val
}


#' Two-sided testing
#'
#' Compute a two-sided p-value given observed test statistic \code{tobs}
#' and randomized \code{tvals}.
#'
#' @inheritParams one_sided_test
#'
#' @details The test may randomize the p-value depending on how many \code{tvals}
#'  are the same as \code{tobs} (up to a tolerance level) to get better performance.
#'
#' @return p-value for two-sided testing.
#' @seealso Testing Statistical Hypotheses (Ch. 15, Lehman and Romano, 2006)
#' @export
two_sided_test = function(tobs, tvals, tol=1e-14) {
  p_val_onesided = one_sided_test(tobs, tvals)
  p_val = 2 * min(p_val_onesided, 1-p_val_onesided)
  p_val
}


#' Generate test statistics
#'
#' Generate default test statistics given a network structure.
#'
#' @param G The \eqn{N\times N} adjacency matrix where \eqn{N} is the number of units. \code{G[i,j]} equals 1
#' iff there is an edge between units \code{i} and \code{j} in the network, and we set
#' all diagonal entries of \code{G} to be 0.
#' @param type The statistics to use. Currently support \code{"elc","score","htn"}
#' in Athey et al. (2018) Exact p-Values for Network Interference section 5.
#' @details In the following \eqn{G} is the adjacency matrix and \eqn{F} is the vector of focal indicator.
#' \itemize{
#' \item \code{"elc"} is the edge-level-contrast test statistic, which equals
#' \eqn{\frac{\sum_{i,j\neq i} F_i G_{ij} (1-F_j) Z_j Y_i^{obs}}{\sum_{i,j\neq i} F_i G_{ij} (1-F_j) Z_j} -
#' \frac{\sum_{i,j\neq i} F_i G_{ij} (1-F_j) (1-Z_j) Y_i^{obs}}{\sum_{i,j\neq i} F_i G_{ij} (1-F_j) (1-Z_j)
#' }}
#'
#' \item \code{"score"} is a score statistic, which equals
#' \eqn{ \mathrm{cov}\left( Y_i^{obs}-\hat{\alpha}-\hat{\tau}_{d} Z_i,~ \sum_{j=1}^N Z_j \bar{G}_{ij} \bigg| \sum_{j=1}^N G_{ij} >0, F_i=1 \right)}.
#'
#' \item \code{"htn"} is the has-treated-neighbor statistic, which equals
#' \eqn{\frac{1}{\mathrm{sd}_{Y_F^{obs}}\mathrm{sd}_{\mathrm{indicator}_F} } \frac{1}{N_F} \sum_{i:F_i=1}
#' (Y_i^{obs}-\bar{Y}_F^{obs})1\{\sum_j G_{ij} Z_j (1-F_j) >0 \} }
#' where \eqn{N_F} is the number of focal units, \eqn{\cdot_F} means restricting a variable to focal units only,
#' and \eqn{\bar{Y}_F^{obs}} is the sample average of \eqn{Y_F^{obs}}.
#' }
#' @export
gen_tstat = function(G, type) {
  stopifnot("type should be one of 'elc', 'score' and 'htn'" = (type %in% c("elc","score","htn")))
  if (type == "elc") {
    f = function(y, z, is_focal) {
      v1 = t(is_focal*y) %*% G %*% ((1-is_focal)*z) / t(is_focal) %*% G %*% ((1-is_focal)*z)
      v2 = t(is_focal*y) %*% G %*% ((1-is_focal)*(1-z)) / t(is_focal) %*% G %*% ((1-is_focal)*(1-z))
      as.numeric(abs(v1-v2))
    }
  }
  if (type == "score") {
    f = function(y, z, is_focal) {
      R = rowSums(G)
      Gij_bar = G/rowSums(G)
      Gij_bar[is.na(Gij_bar)] = 0
      Peer_i = Gij_bar %*% z
      fit0 = lm(y ~ z) # regression model under no peer effects.
      e_r = fit0$residuals
      ### T_score?
      Icond = as.logical(is_focal*(R>0))
      if (!(sum(Icond) > 0)) { return(NaN) }
      as.numeric(abs( cov(e_r[Icond], Peer_i[Icond]) ))
    }
  }
  if (type == "htn") {
    f = function(y, z, is_focal){
      if (sum(is_focal) == 1) { return(NaN) }
      Y_F = y[is_focal]
      sd_F = sd(Y_F)
      neighbor_treated = G %*% (z*(!is_focal)) > 0
      neighbor_treated_F = neighbor_treated[is_focal]
      sd_neighbor_treated_F = sd(neighbor_treated_F)
      Tout = 0
      if ((sd_F>0) & (sd_neighbor_treated_F>0)) {
        Tout = mean( (Y_F-mean(Y_F)) * (neighbor_treated_F-mean(neighbor_treated_F)) ) / (sd_F*sd_neighbor_treated_F)
      }
      as.numeric(abs(Tout))
    }
  }
  return(f)
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


#' Calculating average treatment effect
#'
#' Returns difference in means between exposures \code{b} and \code{a} (as \code{b} - \code{a}), coded as \code{1} and \code{-1}, respectively.
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
  mean(Y[ind1]) - mean(Y[ind2])
}


#' Confidence interval by grid method
#'
#' Compute the CI by inverting the test using a grid method
#'
#' @param ci_lb Pre-determined lower bound of the confidence interval.
#' @param ci_ub Pre-determined upper bound of the confidence interval.
#' @param ci_dec Decimal of the confidence interval.
#'
#' @export
CI_grid = function(ci_lb, ci_ub, ci_dec){
  cat("\n")
  cat("find the confidence interval for the clique-based randomization test using grid method... \n")
  ci_grid = seq(ci_lb, ci_ub, by=10^(-ci_dec)) # generate grids for CI
  #ci = rep(1, length(ci_grid))

  ci = foreach::foreach(igrid = 1:length(ci_grid), .combine='c', .inorder=TRUE) %dopar% {
    tau_temp = ci_grid[igrid]
    Y_temp = out_Yobs(Z_a, Z_b, Y, Y+tau_temp)
    Y_clique_temp = Y_temp[focal_units]
    tobs_temp = ate(conditional_clique[, Zobs_cliqid], Y_clique_temp)
    tvals_temp = apply(as.matrix(conditional_clique[, -Zobs_cliqid]), 2, function(z) ate(z, Y_clique_temp))
    rtest_out_temp = list(tobs=tobs_temp, tvals=tvals_temp)
    decision_temp = out_pval(rtest_out_temp, FALSE, alpha)
    decision_temp
  }

  ci_lb_out = ci_grid[min(which(ci==0))]
  ci_ub_out = ci_grid[max(which(ci==0))]

  return(c(ci_lb_out, ci_ub_out))
}

#' Confidence interval by bisection method
#'
#' Compute the CI by inverting the test using a grid method
#'
#' @inheritParams CI_grid
#'
#' @export
CI_bisection = function(ci_lb, ci_ub, ci_dec){
  cat("\n")
  cat("find the confidence interval for the clique-based randomization test using bisection method... \n")

  # pre-selection of ub
  ci_grid = seq(ci_lb, ci_ub, by=10^(-ci_dec+1))

  stop_ub_init = FALSE
  ub_idx = length(ci_grid)
  ub_init = ci_ub
  while (!stop_ub_init){
    ub_idx = ub_idx - 1
    ub_init_prev = ub_init
    ub_init = tau_temp = ci_grid[ub_idx]
    Y_temp = out_Yobs(Z_a, Z_b, Y, Y+tau_temp)
    Y_clique_temp = Y_temp[focal_units]
    tobs_temp = ate(conditional_clique[, Zobs_cliqid], Y_clique_temp)
    tvals_temp = apply(as.matrix(conditional_clique[, -Zobs_cliqid]), 2, function(z) ate(z, Y_clique_temp))
    rtest_out_temp = list(tobs=tobs_temp, tvals=tvals_temp)
    decision_temp = out_pval(rtest_out_temp, FALSE, alpha)
    if ((decision_temp == 0) | (ub_idx == 1)){stop_ub_init = TRUE}
  }

  # find upper bound
  err = 2 * 10^(-ci_dec)
  ci_ub_bis = ub_init_prev
  ci_lb_bis = ub_init
  ci_ub_out = (ci_ub_bis + ci_lb_bis) / 2
  while (err>10^(-ci_dec-1)){ # err is one order less than ci_dec, so that we can round it by ci_dec
    tau_temp = ci_ub_out
    Y_temp = out_Yobs(Z_a, Z_b, Y, Y+tau_temp)
    Y_clique_temp = Y_temp[focal_units]
    tobs_temp = ate(conditional_clique[, Zobs_cliqid], Y_clique_temp)
    tvals_temp = apply(as.matrix(conditional_clique[, -Zobs_cliqid]), 2, function(z) ate(z, Y_clique_temp))
    rtest_out_temp = list(tobs=tobs_temp, tvals=tvals_temp)
    decision_temp = out_pval(rtest_out_temp, FALSE, alpha)
    if (decision_temp == 0){
      ci_ub_bis = ci_ub_bis
      ci_lb_bis = ci_ub_out
    } else {
      ci_ub_bis = ci_ub_out
      ci_lb_bis = ci_lb_bis
    }
    err = abs(ci_ub_out - (ci_ub_bis + ci_lb_bis) / 2)
    ci_ub_out = (ci_ub_bis + ci_lb_bis) / 2
  }

  # pre-selection of lb
  stop_lb_init = FALSE
  lb_idx = 1
  lb_init = ci_lb
  while (!stop_lb_init){
    lb_idx = lb_idx + 1
    lb_init_prev = lb_init
    lb_init = tau_temp = ci_grid[lb_idx]
    Y_temp = out_Yobs(Z_a, Z_b, Y, Y+tau_temp)
    Y_clique_temp = Y_temp[focal_units]
    tobs_temp = ate(conditional_clique[, Zobs_cliqid], Y_clique_temp)
    tvals_temp = apply(as.matrix(conditional_clique[, -Zobs_cliqid]), 2, function(z) ate(z, Y_clique_temp))
    rtest_out_temp = list(tobs=tobs_temp, tvals=tvals_temp)
    decision_temp = out_pval(rtest_out_temp, FALSE, alpha)
    if ((decision_temp == 0) | (lb_idx == length(ci_grid))){stop_lb_init = TRUE}
  }

  # find lower bound
  err = 2 * 10^(-ci_dec)
  ci_ub_bis = lb_init
  ci_lb_bis = lb_init_prev
  ci_lb_out = (ci_ub_bis + ci_lb_bis) / 2
  while (err>10^(-ci_dec-1)){
    tau_temp = ci_lb_out
    Y_temp = out_Yobs(Z_a, Z_b, Y, Y+tau_temp)
    Y_clique_temp = Y_temp[focal_units]
    tobs_temp = ate(conditional_clique[, Zobs_cliqid], Y_clique_temp)
    tvals_temp = apply(as.matrix(conditional_clique[, -Zobs_cliqid]), 2, function(z) ate(z, Y_clique_temp))
    rtest_out_temp = list(tobs=tobs_temp, tvals=tvals_temp)
    decision_temp = out_pval(rtest_out_temp, FALSE, alpha)
    if (decision_temp == 0){
      ci_ub_bis = ci_lb_out
      ci_lb_bis = ci_lb_bis
    } else {
      ci_ub_bis = ci_ub_bis
      ci_lb_bis = ci_lb_out
    }
    err = abs(ci_lb_out - (ci_ub_bis + ci_lb_bis) / 2)
    ci_lb_out = (ci_ub_bis + ci_lb_bis) / 2
  }

  if (ub_idx == 1){
    cat("The upper bound is below the pre-set lower bound")
    ci_ub_out = ci_lb
  } else if (ub_idx == length(ci_grid) - 1){
    cat("The upper bound is above the pre-set upper bound")
    ci_ub_out = ci_ub
  } else {ci_ub_out = round(ci_ub_out, ci_dec)} # err is one order less than ci_dec, hence can round it by ci_dec

  if (lb_idx == 2){
    cat("The lower bound is below the pre-set lower bound")
    ci_lb_out = ci_lb
  } else if (lb_idx == length(ci_grid)){
    cat("The lower bound is above the pre-set upper bound")
    ci_lb_out = ci_ub
  } else {ci_lb_out = round(ci_lb_out, ci_dec)}

  return(c(ci_lb_out, ci_ub_out))
}

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

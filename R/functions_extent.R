#' Test function for "extent of interference" null hypotheses.
#'
#' @param Y The observed outcome vector.
#' @param Z A binary matrix of dimension (number of units x number of randomizations, i.e. assignments.) storing the assignment vectors.
#' @param expos A three-way array with dimensions (number of units x number of randomizations x dimension of exposure).
#' Its entry records the exposure \code{f_i(z)} of units under treatments, where the exposure has more than one dimensions.
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}.
#' @param decom The algorithm used to calculate the biclique decomposition. Currently supported algorithms are "bimax" only.
#' @param ... Other stuff ...
#'
#' @return A list of items summarizing the randomization test.
clique_test_ex = function(Y, Z, expos, Zobs_id, decom, ...){

  addparam = list(...) # catch variable parameters

  # setting default values if they do not exist
  if(is.null(addparam$exclude_treated)){ exclude_treated=TRUE } else { exclude_treated=addparam$exclude_treated }
  if(is.null(addparam$stop_at_Zobs)){ stop_at_Zobs=FALSE } else { stop_at_Zobs= addparam$stop_at_Zobs}
  if(is.null(addparam$ret_pval)){ ret_pval=TRUE } else { ret_pval=addparam$ret_pval }
  if(is.null(addparam$adj_Y)){ adj_Y=FALSE } else { adj_Y=(addparam$adj_Y)&(!is.null(Xadj)) } # is TRUE iff we specify it and Xadj is not null

  # decompose the null-exposure graph
  cat("decompose the null-exposure graph ... \n")
  Z0 = 1:dim(Z)[2]
  conditional_clique_idx = list(focal_unit=c(), focal_ass=c())
  conditional_clique = Matrix(0, nrow=dim(Z)[1], ncol=dim(Z)[2], sparse=TRUE)
  if (decom == 'bimax'){
    while (length(Z0)>0){
      cat("\r","Z0 now has length", length(Z0), '...')
      # 1. select one random from 1:length(Z0) as Zsub
      # 2. generate multiNEgraph = out_NEgraph_ex(Zsub, Z0, expos)
      # 3. decompose multiNEgraph to multi_clique, get one biclique is enough
      # 4. conditional_clique =union of multi_clique, delete multi_clique's col from Z0

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

  # greedy decomposition is still under development
  # if (decom == 'greedy'){
  #   decomp = out_clique_decomposition_greedy(NEgraph, Zobs_id, num_ass=addparam$minass, stop_at_Zobs)
  #   if(stop_at_Zobs){
  #     conditional_clique = decomp
  #   } else {conditional_clique = out_clique(Zobs_id,decomp)}
  # }

  # test
  cat("\n finding test statistics ... \n")
  conditional_clique = as.matrix(conditional_clique)
  teststat = apply(conditional_clique, 2, function(z) {mean(Y[as.logical(z)])})
  tobs = teststat[Zobs_id]; tvals = teststat[setdiff(conditional_clique_idx$focal_ass, Zobs_id)]
  pval = out_pval(list(tobs=tobs, tvals=tvals), T, 0.05) # here we use the previous out_pval function

  retlist = list(pval=pval, tobs=tobs, tvals=tvals, conditional_clique=conditional_clique)
  return(retlist)

}


#' Generate the multi-null exposure graph.
#'
#' @param Zsub An integer from \code{1:length(Z0)} that represents a sample from \code{Z0}.
#' @param Z0 A vector of column index that represents the remaining treatment.
#' Please see Sec. 8 of the paper for more details.
#' @param expos A three-way array with dimensions (number of units x number of randomizations x dimension of exposure).
#' Its entry records the exposure \code{f_i(z)} of units under treatments, where the exposure has more than one dimensions.
#'
#' @return A matrix of dimension (number of units x \code{length(Z0)}) which is the multi-null exposure graph with
#' respect to \code{Z0[Zsub]} and \code{Zsub}.
out_NEgraph_ex = function(Zsub, Z0, expos, tol=1.5e-8){
  expos_sub = expos[, Z0, ]
  num_unit = dim(expos_sub)[1]
  NEgraph = Matrix(0, nrow=num_unit, ncol=length(Z0), sparse=TRUE)

  # subtract each column from the Zsub column
  expos_sub_idx = apply(expos_sub, 3, function(expos_1d) {abs(sweep(expos_1d, 1, expos_1d[,Zsub])) < tol}, simplify = F)
  expos_sub_idx = Reduce("*", expos_sub_idx)
  NEgraph[expos_sub_idx == 1] = 1

  rownames(NEgraph) = 1:num_unit
  colnames(NEgraph) = 1:length(Z0)
  return(as.matrix(NEgraph))
}




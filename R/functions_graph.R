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
  expos_sub_Zsub = expos_sub[,Zsub,]
  num_unit = dim(expos_sub)[1]; num_rand = dim(expos_sub)[2]; dim_expos_sub = dim(expos_sub)[3]

  # subtract each column from the Zsub column
  for (id in 1:num_rand){
    expos_sub[,id,] = expos_sub[,id,] - expos_sub_Zsub
  }
  NEgraph = (rowSums(abs(expos_sub)<tol, dims=2) == dim_expos_sub) * 1

  rownames(NEgraph) = 1:num_unit
  colnames(NEgraph) = 1:length(Z0)
  return(NEgraph)
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
#' @param decomp A list containing the clique decomposition of the null-exposure graph, given by one of the decomposition algorithms.
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


#' Decomposing Null Exposure Graph Using Bimax
#'
#' One of the main functions for implementing the methodology.
#' Outputs a biclique decomposition of the null-exposure graph using Bimax algorithm.
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
out_clique_decomposition_bimax = function(NEgraph,Zobs_id,minr,minc,stop_at_Zobs=TRUE){

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

#' Greedy Algorithm to Find a Biclique
#'
#' This function returns a biclique in a null-exposure graph using a greedy algorithm. The
#' minimum number of focal assignments in the biclique is \code{num_ass}.
#'
#' @param NEgraph The null-exposure graph we want to decompose.
#' @param num_ass Controls the minimum number of focal assignments in the biclique found
#' by the algorithm (algorithm runtime is sensitive to this value).
#' @return A biclique of the input null-exposure graph \code{NEgraph}.
#' @export
out_greedy_decom = function(NEgraph, num_ass){
  focal_units = c() # set of focal units
  focal_ass = c() # set of focal assignments

  CONT = TRUE
  NEgraph_original = NEgraph
  get_clique_size = function(units) {
    ne_temp = matrix(NEgraph_original[units,], ncol=ncol(NEgraph_original))
    # sum(apply(ne_temp^2, 2, prod)) # the apply() gives 1 for assignments which all focals are connected to.
    sum(colSums(ne_temp) == length(units))
  }
  get_clique = function(units) {
    ne_temp = matrix(NEgraph_original[units,], ncol=ncol(NEgraph_original))
    # a = apply(ne_temp^2, 2, prod) # =1 for assignments which all focals are connected to.
    a = (colSums(ne_temp) == length(units)) * 1 # should not use as.numeric: keep the column names
    keep_ass = which(a==1)
    return(NEgraph_original[units, keep_ass])
  }

  # if the NEgraph contains a biclique that incorporates all units, and numcol>num_ass, return directly
  if (get_clique_size(1:dim(NEgraph_original)[1]) > num_ass ) {
    return(list(clique=get_clique(1:dim(NEgraph_original)[1])))
  }

  while(CONT) {
    # print(paste("Dim of NEgraph=",dim(NEgraph)))
    units = as.numeric(rownames(NEgraph))  # 2, 5, 6, ...
    ass = as.numeric(colnames(NEgraph))

    # degree = as.numeric(rowSums(NEgraph))
    # rand_id = sample(1:length(units), 1, prob = degree) # random unit
    rand_id = sample(1:length(units), 1) # select a random unit, only an index of that unit but not rowname.
    focal_unit = units[rand_id]  # rowname of that randomly selected unit.

    # print(paste("Select", focal_unit)) # unit name.
    iass_remove = which(NEgraph[rand_id,]==0)  # remove assignments that are not connected to this randomly selected unit.
    NEgraph = NEgraph[-c(rand_id),-iass_remove]
    focal_units = c(focal_units, focal_unit)  # add this unit to the set of focal units.

    nsize = get_clique_size(focal_units)
    if(nsize > num_ass) {
      # print(paste("Added focal. Total", length(focal_units),"focal units. Clique size=", nsize,'\n'))
    } else { # the algorithm ends if after adding this unit, num of assignments in the biclique is below num_ass
      focal_units = setdiff(focal_units, focal_unit)
      CONT = FALSE
    }
  }

  focal_ass = as.numeric(colnames(NEgraph))
  clique = get_clique(focal_units)
  return(list(clique=clique))
}

#' Decomposing Null Exposure Graph Using the Greedy Algorithm
#'
#' Outputs a biclique decomposition of the null-exposure graph using the greedy algorithm
#' proposed above. The resulting decomposition paritions the assignment space,
#' while the unit space may overlap.
#'
#' @param NEgraph The null-exposure graph
#' @param Zobs_id The index location of the observed assignment vector in \code{Z}, \code{Z_a}, and \code{Z_b}.
#' @param num_ass Controls the minimum number of focal assignments in the biclique found
#' by the algorithm (algorithm runtime is sensitive to this value).
#' @param stop_at_Zobs A Boolean indicating whether the decomposition algorithm should stop when the clique containing the observed assignment is found.
#' @return If \code{stop_at_Zobs} is \code{TRUE}, a matrix representing the clique to condition upon.
#' If \code{stop_at_Zobs} is \code{FALSE}, a list containing the clique decomposition of the null-exposure graph.
#' @export
out_clique_decomposition_greedy = function(NEgraph, Zobs_id, num_ass, stop_at_Zobs){

  iremove = which(rowSums(NEgraph!=0)==0)  # removes isolated units.
  if(length(iremove)!=0){ NEgraph = NEgraph[-iremove,] }

  cc = 1
  CONT = TRUE
  allcliques = list()

  while(CONT){
    test = out_greedy_decom(NEgraph, num_ass)
    clique = test$clique
    focalass = colnames(clique)
    zobs_there = which(focalass==Zobs_id)
    if(stop_at_Zobs){ if(length(zobs_there)!=0){ CONT=FALSE } } # see whether selected clique contains Zobs
    if(!stop_at_Zobs){  }
    iass_remove = match(focalass,colnames(NEgraph))
    NEgraph = NEgraph[,-iass_remove]
    cat('\r','found greedy clique',cc,'...')
    allcliques[[cc]] = clique
    if(dim(NEgraph)[2]<=num_ass){ # when remaining cols not enough to do the greedy algo.
      units_leftover = which(rowSums(NEgraph^2)==dim(NEgraph^2)[2])
      allcliques[[cc+1]] = NEgraph[units_leftover,]
      CONT=FALSE
    }
    cc = cc+1
  }

  if(stop_at_Zobs){
    returnlist = clique
  } else {returnlist = allcliques}

  return(returnlist)
}

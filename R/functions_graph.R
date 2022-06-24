#' Generate decompositions of multi-null exposure graphs given null hypothesis
#'
#' @inheritParams clique_test
#'
#' @return a list \code{MNE} of clique decomposition and specified \code{controls}.
#' Each element of \code{MNE} records one biclique decomposed from a multi-null exposure
#' graph and contains its focal units and focal assignments.
#'
#' @export
biclique.decompose = function(Zrealized, hypothesis,
                              controls=list(method="greedy", mina=10, num_randomizations=2000)){

  # catch functions and controls from hypothesis and controls
  design_fn = hypothesis$design_fn
  exposure_i = hypothesis$exposure_i
  null_equiv = hypothesis$null_equiv
  num_randomizations = controls$num_randomizations
  if (is.null(num_randomizations)) {num_randomizations = 2000} # set to be 2000 if not supplied

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
    mina = controls$mina
    if (!is.numeric(mina)) {
      stop("if 'greedy' is used in controls, should supply 'mina' as integer", call. = F)
    }
  } else {
    stop("the decomposition method in controls should be either 'bimax' or 'greedy'", call. = F)
  }

  # generate randomizations using design_fn & num_randomizations --> Z_m
  num_units = length(Zrealized)
  Z_m = matrix(0, nrow=num_units, ncol=(num_randomizations+1))
  Z_m[, 1] = Zrealized
  for (id_rand in 1:num_randomizations){
    Z_m[, id_rand+1] = design_fn()
  }
  Zobs_id = 1

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
        num_ass = mina
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

  return(list(MNE=MNE, controls=controls, method=method))
}



#' Generate the multi-null exposure graph.
#'
#' @param Zsub An integer from \code{1:length(Z0)} that represents a sample from \code{Z0}.
#' @param Z0 A vector of column index that represents the remaining treatment.
#' Please see Sec. 8 of the paper for more details.
#' @param Z A matrix of dimension (\code{num_units} x \code{num_randomizations+1}) that records the randomizations
#' of treatments for each unit.
#' @param num_units Number of units in the experiment, equals to the length of the outcome vector.
#' @param exposure_fn A function that specifies the equality of exposure
#' of a pair of individuals under different randomization.
#'
#' @return A matrix of dimension (\code{num_units} x \code{length(Z0)}) which is the multi-null exposure graph with
#' respect to \code{Z0[Zsub]} and \code{Z0}.
out_NEgraph_multi = function(Zsub, Z0, Z, num_units, exposure_fn){
  multiNEgraph = matrix(0, nrow = num_units, ncol = length(Z0))
  Zsub_treatment = Z[,Z0[Zsub]]
  for (i in 1:num_units){
    for (zid in 1:length(Z0)){
      z = Z[,Z0[zid]]
      multiNEgraph[i, zid] = exposure_fn(Zsub_treatment, z, i)
    }
  }

  multiNEgraph = multiNEgraph * 1
  rownames(multiNEgraph) = 1:num_units
  colnames(multiNEgraph) = 1:length(Z0)
  return(multiNEgraph)
}

#' Generate multi-null exposure graph using null_equiv.
#'
#' @inheritParams clique_test
#' @param Zsub An integer from \code{1:length(Z0)} that represents a sample from \code{Z0}.
#' @param Z0 A vector of column index that represents the remaining treatment.
#' Please see Sec. 8 of the paper for more details.
#' @param expos A three-way array with dimensions (number of units x number of randomizations x dimension of exposure).
#' Its entry records the exposure \code{f_i(z)} of units under treatments, where the exposure has more than one dimensions.
#'
#' @return A matrix of dimension (\code{num_units} x \code{length(Z0)}) which is the multi-null exposure graph
#' with respect to \code{Z0[Zsub]} and \code{Z0}.
out_NEgraph_multi_separate = function(Zsub, Z0, expos, null_equiv, dim_exposure){
  expos_sub = expos[, Z0, ]
  num_units = dim(expos_sub)[1]; num_rand = dim(expos_sub)[2]; dim_expos_sub = dim(expos_sub)[3]

  multiNEgraph = matrix(0, nrow = num_units, ncol = num_rand)
  if (dim_exposure == 1){
    for (i in 1:num_units){
      # compare equivalence of Zsub and each z in Z0
      for (zid in 1:num_rand){
        multiNEgraph[i, zid] = null_equiv(expos_sub[i, zid, 1], expos_sub[i, Zsub, 1])
      }
    }
  } else {
    for (i in 1:num_units){
      # compare equivalence of Zsub and each z in Z0
      for (zid in 1:num_rand){
        multiNEgraph[i, zid] = null_equiv(expos_sub[i, zid, ], expos_sub[i, Zsub, ])
      }
    }
  }

  multiNEgraph = multiNEgraph * 1
  rownames(multiNEgraph) = 1:num_units
  colnames(multiNEgraph) = 1:length(Z0)
  return(multiNEgraph)
}


#' Generating the Null Exposure Graph for contrast hypothesis.
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
out_NEgraph_contrast = function(Z_a,Z_b,Z,exclude_treated=TRUE){
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
#' @param NEgraph The null-exposure graph.
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
    # the input units are actual row names of NEgraph_original
    units
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

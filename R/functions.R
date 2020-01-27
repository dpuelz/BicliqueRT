#' Load in necessary libraries
get_libraries = function(){
  library(readstata13)
  # library(biclique)
  library(Matrix)
  library(filematrix)
  library(data.table)
  library(pbapply)
  library(lattice)
  library(grid)
  library(gridExtra)
}

#' compute the observed outcome vector, given potential outcomes on the exposures
get_Yobs = function(Z,ya,yb){
  Za = (Z==-1)
  Zb = (Z==1)
  y = ya*Za + yb*Zb
  y
}

get_clique = function(z.id,decom){
  for(ii in 1:length(decom)){
    if(sum(z.id==colnames(decom[[ii]]))){
      return(ii)
    }
  }
}

get_bicliquedecom = function(NEmat,Zobs_id,minr=50,minc=50,stopatZobs=FALSE){
  library(biclust)

  iremove = which(rowSums(NEmat!=0)==0)  # removes isolated units.
  if(length(iremove)!=0){ NEmat = NEmat[-iremove,] }
  numb = 1
  ncol = ncol(NEmat)

  decomp = c()
  new.NEmat = NEmat[,1:ncol]
  numleft = ncol(new.NEmat)
  oldnames = colnames(new.NEmat)
  ii=1

  while(numleft>0){
    minc.new = min(minc,numleft)
    bitest = biclust(new.NEmat, method=BCBimax(), minr, minc.new, number=numb)
    bicliqMat = bicluster(new.NEmat,bitest)
    themat = bicliqMat$Bicluster1
    while(length(themat)==0){
      numleft=numleft-1
      minc.new = min(minc,numleft)
      bitest = biclust(new.NEmat!=0, method=BCBimax(), minr, minc.new, number=numb)
      bicliqMat = bicluster(new.NEmat,bitest)
      themat = bicliqMat$Bicluster1
    }
    dropnames = colnames(themat)
    oldnames = colnames(new.NEmat)
    drop.ind = match(dropnames,oldnames)

    if(stopatZobs){
      if(sum(dropnames==Zobs_id)==1){
        themat.Zobs.s = themat; cat("found clique with Zobs!\n")
        return(themat.Zobs.s)
      }
    }
    new.NEmat = new.NEmat[,-drop.ind]
    numleft = dim(new.NEmat)[2]
    if(length(numleft)==0){ numleft=0 }
    decomp[[ii]] = themat
    cat("found clique",ii,'...\n')
    ii=ii+1
  }
  return(decomp)
}

run_test = function(Yobs,BImat,Zobs_id){
  focalunits = as.numeric(rownames(BImat))
  cliqassign = as.numeric(colnames(BImat))

  # restrict Y to clique
  Y.clique = Yobs[focalunits]

  # test statistics
  Zobs_cliqid = which(Zobs_id==cliqassign)
  Tobs = ate(BImat[, Zobs_cliqid], Y.clique)

  Trand = apply(as.matrix(BImat[, -Zobs_cliqid]), 2, function(z) ate(z, Y.clique))

  # pval = mean(Trand>=Tobs)
  reject = (one_sided_test(tobs=Tobs,tvals=Trand,alpha=0.05,tol=1e-14))
  reject
}

create_network = function(numnodes){
  # generate new network (3 clusters of 2D Gaussians)
  x = c(rnorm(numnodes*0.5,sd=0.1)+0.5,rnorm(numnodes*0.3,sd=0.075)+0.25,rnorm(numnodes*0.2,sd=0.075)+0.3)
  y = c(rnorm(250,sd=0.1)+0.5,rnorm(150,sd=0.075)+0.75,rnorm(100,sd=0.075)+0.3)
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

  save(Dmat,file=paste('Dsim_nodes=',numnodes,'.RData',sep=''))
  thepoints
}

sparsify <- function(mat){
  mat = Matrix(mat,sparse=TRUE)
  mat
}

get_DZ2 = function(r,D,Z){
  Dr = binarize_r2(r,D)
  Z = sparsify(Z)
  DZ = Dr%*%Z
  DZ
}

get_binaryDZ = function(DZ,direction){
  if(direction=='lt'){
    DZbin = (DZ>0)
  }
  if(direction=='gt'){
    DZbin = (DZ==0)
  }
  DZbin
}

binarize_r2 <- function(r,mat){
  # r is a distance metric, such as euclidean distance between street segments
  # every segment pair > r distance apartment will be zero.
  matnew = Matrix(0,dim(mat)[1],dim(mat)[2],sparse=TRUE)
  rownames(matnew) = rownames(mat)
  colnames(matnew) = colnames(mat)

  matnew[mat<=r] <- 1
  diag(matnew) <- -1e6 # this is key here .. focus hypothesis on untreated units only
  matnew
}

ate = function(Z,Y){
  ind1 = which(Z==1)
  ind2 = which(Z==-1)
  mean(Y[ind1])-mean(Y[ind2])
}

one_sided_test = function(tobs, tvals, alpha, tol=1e-14) {
  tvals = c(tobs,tvals)
  srt = sort(tvals)
  M = length(tvals)
  k = ceiling(M * (1-alpha))
  Tk = srt[k]
  if(abs(tobs - Tk) < tol) {
    # if tobs=Tk
    ax = (M * alpha - sum(tvals > Tk)) / sum(abs(tvals - Tk) < tol)
    return(runif(1) <= ax) ## randomize.
  }

  return(tobs > Tk)
}

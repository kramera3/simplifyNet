# Author: Andrew M. Kramer
# Edited: Alexander M. Mercier [2020]

require(igraph)
require(cPCG)
require(Matrix)

source('Matrix_Mod.R')
source('SpielmanSrivastava_Sparse.R')


# Script for network sparsification
# Will have a single wrapper function that calls options for each version
# Will rely on supplied functions from outside for refitting, returning negative log-likelihood (or similar), these functions should be sourced and available in the environment
# User also supplies options

# Network Sparsification
# Possible methods
  # 1. Toivonen - Remove all edges that are not part of best path, algorithm from Toivonen et al. implemented with 'igraph' package
  # 2. GlobalSparse - Remove edges under an edge weight cutoff
  # 3. LocalAdapt - Remove weakest edges locally rather than globally, can set proportion or use refitting and cutoff
  # 4. EffectiveResistance - Randomly samples edges with probability proportional to the effective resistance of that edge, algorithm from Spielman and Srivastava (untested on directed networks)

simplifyNet <- function(data, method="Toivonen", model, func, cutoff, remove.prop, num.samples, epsilon, matrix.sparse = FALSE, num.nodes = NULL){
  # Removes edges from network models using algorithm of choice
  #
  # Args:
  #   data: edge list structured | node 1 | node 2 | weight | as class data.frame, adjacency matrix of class matrix or sparseMatrix object
  #   model: model to run on network for iterative sparsification process (NOT YET IMPLEMENTED)
  #   func: scoring function for iterative sparsification (NOT YET IMPLEMENTED)
  #   method: algorithm used to remove edges, currently includes Toivonen (best-path), Refit and LocalAdapt. Default is Toivonen.
  #   cutoff: minimum edge weight to retain in simplified network
  #   search: whether to remove chunks of edges or proceed one-at-a-time. "direct" refits after each edge, can often be computationally prohibitive
  #   disconnect: whether to allow sparsification to lead to disconnection, only relevant for GlobalSparse
  #   remove.prop: argument can be passed to LocalAdapt or GlobalSparse to remove a set proportion of edges instead of using statistical cutoff
  #   num.samples: number of samples to take for Effective Resistance sampling
  #   epsilon: level of approximation for calculating effective resistance
  #   (to ensure quadratic laplacian similarity to multiplicative error alpha, let q = 9*C^2*n*logn/alpha^2 where C is a suitably large absolute constant for concentration)
  #   matrix.sparse: choose to return sparseMatrix adjacency matrix
  #   num.nodes: if there are disconnected nodes, specify the number of nodes for output
  # Output:
  #   output: sparsifed network in the class as inputted data

  #Ensure that data is in edge list format if input data is an adjacency matrix
  sm.lock <- 0
  mtrx.lock <- 0
  if (!is.null(num.nodes)  & !is.matrix(data)){
    num.nodes <- max(data[,1:2])
  } else {
    num.nodes <- dim(data)[1]
  }

  if (is(data, 'sparseMatrix')) {
    sm.lock <- 1
    data <- data.frame(Mtrx_Elist(as.matrix(data)))
  }
  if (method == "Toivonen"){
    if (is.data.frame(data) == T){
      data < - Elist_Mtrx(data, num.nodes)
    }
    mtrx.lock <- 1
    sparse <- toivonen(data=data)
  }
  if (method == "GlobalSparse"){
    if (is.matrix(data) == T){
      data <- data.frame(Mtrx_Elist(data))
      mtrx.lock <- 1
    }
    sparse <- gns(data=data, cutoff=cutoff, num.nodes=num.nodes)
  }
  if (method == "LocalAdapt"){
    if (is.matrix(data) != T){
      data <- Elist_Mtrx(data, num.nodes)
    } else {
      mtrx.lock <- 1
    }
    sparse <- lans(data=data, remove.prop=remove.prop)
  }
  if (method == 'EffectiveResistance'){
    if (is.matrix(data) != T){
      data <- Elist_Mtrx(data, num.nodes)
    } else {
      mtrx.lock <- 1
    }
    R <- EffR(network = data, epsilon = epsilon, n = NULL) #High RAM usage...
    sparse <- EffRSparse(network = data, q = num.samples, R = R)
  }

  #Class of value inputted is outputted
  if (isTRUE(matrix.sparse) | sm.lock == 1) {
    output <- as(sparse, 'sparseMatrix')
  } else if (mtrx.lock == 1){
    output <- sparse
  } else {
    output <- data.frame(Mtrx_Elist(sparse))
    colnames(output) <- c("node 1", "node 2", "weight")
  }

  return(output)
}

# Returns the sparse network containing every link that is in a best path, Toivonen et al
toivonen <- function(data){
  ## Determine if network is symmetric or not ##
  if (isSymmetric.matrix(data)){
    mode = 'undirected'
  } else {
    mode = 'directed'
  }

  grph <- graph_from_adjacency_matrix(data, mode = mode, weighted = TRUE)
  E(grph)$weight <- (E(grph)$weight)^(-1) #Change edge weights for best path

  n <- dim(data)[1]
  ## Make solution matrix ##
  mask<-matrix(FALSE, nrow=n, ncol=n)
  ## Look for best paths (in each direction) ##
  for (i in 1:n){
    paths <- get.shortest.paths(grph, from=as.character(i))[[1]] #The initial list has 4 elements, all the best paths are stored in element[[1]] #Does this include weights?
    ## Transform path to "mask" indices of 1 ##
    for(j in 1:length(paths)){
      if(length(paths[[j]]) > 1){ #The path of each to itself is undefined, numeric(0)
        for(k in 1:(length(paths[[j]]) - 1)){
          mask[paths[[j]][k], paths[[j]][k+1]] <- TRUE
        }
      }
    }
  }

  return(mask * data)
}

#Returns the sparse network by some global threshold
gns <- function(data, cutoff, num.nodes){
  if (cutoff > 0) {
    colnames(data)<-c("n1","n2","weight")
    ret <- data[data$weight %in% data$weight[data$weight > cutoff],]
    ret <- Elist_Mtrx(ret, num.nodes)

    return(ret)
  }
}

#Returns the sparse network based on Locally adaptive network sparsification, Foti et al 2011, modified for directed network
lans <- function(data, remove.prop){
  n <- dim(data)[1]
  Aij<-matrix(0, nrow=n, ncol=n)
  ##Get the fractional weights of the edges for each node
  ##go across rows and divide by sum, then across columns, calculating one triangle at a time
  ##In the models source is j and destination is i, that explains switch of rows and columns with Foti et al Supplement S4
  Pij<-prop.table(data,margin=1) #Relative prob of numeric in the matrix across rows.
  lower<-prop.table(data, margin=2) #Relative prob of numeric in the matrix across columns. #Note# This will be the same as above iff the matrix is symmetric
  Pij[lower.tri(Pij)]<-lower[lower.tri(lower)]#lower.tri gives a matrix of same size with TRUE in the lower triangle (default is not to include diagonal)
  if (length(remove.prop>0)){
    for (i in 1:n){
      temp<-ecdf(Pij[,i])
      pa<-quantile(temp,1-remove.prop)
      for (j in 1:n){
        if (Pij[j,i]>pa) Aij[j,i]<-1
      }
    }
  } else { # call the iterative refitting algorithm that checks against cutoff
    Aij<-Pij
  }

  return(Aij*data)
}


toy_iterative<-function(data,remove.prop,remove.per){
  if (is.matrix(data) == TRUE) {
    data <- adj_to_edgelst(data) #Convert to edge list
    colnames(data) <- c("from", "to", "weight")
  }
  it <- temp <- data[order(data$weight),] #Order edge list
  o_w <- sum(data$weights)
  while(TRUE) {
    if(count_components(graph_from_data_frame(temp))>1){ #Check if network is connected
      return(it)
    }else{
      it <- temp
      n_w <- sum(temp$weight) #Compute new sum of weights
      if(n_w < remove.prop * o_w){ #Check if new sum is below remove.prop
        return(it)
      } else {
        n = ceiling(remove.per*length(it)) #Remove percent of edge list
        for(i in 1:n){
          temp <- it[-i,]
        }
      }
    }
  }
}

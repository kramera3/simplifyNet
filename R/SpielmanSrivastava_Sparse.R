#' @author Author: Alexander M. Mercier
#' @rdname SpielmanSrivastava_Sparse

# From Spielman and Srivastava 2011

#source('Matrix_Mod.R')

#' @name Qgen
#' @title Random Bernoulli Matrix Generation
#' @description  Compute random +/-\code{1/sqrt(k)} Bernoulli matrix where \code{k=24log(n)/epsilon^2} \cr
#' Used to fulfill Johnson-Lindenstrauss lemma \cr
#' From Spielman and Srivastava 2011 \cr
#  Input:
#' @param  n
#' Number nodes in the given network
#' @param epsilon
#' Used to specify resolution of the approximation with regard to the Johnson-Lindenstrauss lemma
# Output:
#' @return Generates a random +/- 1/sqrt(k) Bernoulli matrix, \code{Q}
#' @export
QGen <- function(n, epsilon){
  k = ((log(n,2))) / (epsilon)
  k_d = ceiling(k)
  x = 1 / sqrt(k)
  Q = matrix(sample(c(-x,x), size = k_d * n, replace = TRUE),
             nrow = k_d,
             ncol = n)

  return(Q)
}


#' @name EffR
#' @title Effective Resistance Approximation
#' @description  Approximately calculate the effective resistances between all nodes of a graph. \cr
#' Based on algorithm from Spielman-Srivastava (2011) \cr
#' This version of the effective resistance approximation calculator is RAM intensive. If you are working with a large network, use either \link[simplifyNet]{EffR2} or \link[simplifyNet]{EffRPar}
#' @import igraph cPCG
# Input:
#' @param network
#' An adjacency matrix or edge list of the network; edge list should be in the format |in|out|weight|
#' @param epsilon
#' Parameter to control resolution of the approximation
#' @param n
#' The number of nodes; default is \code{NULL} and should not be changed unless the network is disconnected
# Output:
#' @return A matrix, \code{R}, such that \code{R_ij} gives the approximate effective resistance between node \code{i} and \code{j}
#' @examples
#' A <- ER_gen(n=100, p=0.1, weights=1)
#' R <- EffR(network=A, epsilon=0.1)
#' @export
EffR <- function(network, epsilon, n = NULL){
  if (nrow(network) != ncol(network)){
    if (is.null(n)) n = max(network[,1:2])
    network = EList_Mtrx(network, n)
  }
  A = methods::as(network, 'matrix')
  B = sVIM(network)
  L = as.matrix(diag(rowSums(abs(network))) - network)
  W = WDiag(network)
  Q = QGen(nrow(B), epsilon) #How to make this less RAM intensive?
  Y = methods::as((Q %*% (sqrt(W)) %*% B), 'matrix')
  rm(Q, W, B)
  Z = matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
  for (i in 1:nrow(Y)){
    Y_v = matrix(Y[i, ])
    Z[i,] = cPCG::cgsolve(L, Y_v)
  }
  rm(Y)
  R = matrix(0L, nrow = nrow(L), ncol = ncol(L))
  for (i in 1:nrow(L)){
    for (j in 1:nrow(L)){
      if (i != j) R[i, j] = abs(sum((Z[,i] - Z[,j])^2))
    }
  }

  return(R)
}


#' @name EffR2
#' @title Effective Resistance Approximation with RAM Features
#' @description  Approximately calculate the effective resistances between all nodes of a graph. \cr
#' Based on algorithm from Spielman-Srivastava (2011) \cr
#' This version of the effective resistance calculator is more RAM efficient, but slower than \link[simplifyNet]{EffR}; if speed is a concern and parallel processing is possible, \link[simplfyNet]{EffRPar} is likely a better alternative
#' @import igraph cPCG
# Input:
#' @param network
#' An adjacency matrix or edge list of the network; edge list should be in the format |in|out|weight|
#' @param epsilon
#' Parameter to control resolution of the approximation
#' @param n
#' The number of nodes; default is \code{NULL} and should not be changed unless the network is disconnected
# Output:
#' @return A matrix, \code{R}, such that \code{R_ij} gives the approximate effective resistance between node \code{i} and \code{j}
#' @examples
#' A <- ER_gen(n=100, p=0.1, weights=1)
#' R <- EffR(network=A, epsilon=0.1)
#' @export
EffR2 <- function(network, epsilon, n = FALSE){
  if (nrow(network) != ncol(network)){
    if (is.null(n)) n = max(network[,1:2])
    network = EList_Mtrx(network, n)
  } else {
    n = nrow(network)
  }
  m = dim(Mtrx_EList(network))[1]
  tolProb = 0.5
  L = as.matrix(diag(rowSums(abs(network))) - network)
  B = sVIM(as.matrix(network))
  W = WDiag(as.matrix(network))
  k = ceiling((24*log(n,2)) / (epsilon^2))
  Z = matrix(NA, nrow=k, ncol=n)
  for(i in 1:k){
    ons = stats::runif(m, 0, 1) > tolProb
    ons = ons - !(ons)
    ons = ons/sqrt(k)
    Z_r = cgsolve(L,t(as.matrix(ons%*%W%*%B)))
    Z_r = t(Z_r)
    Z[i,] = Z_r
  }
  R = matrix(0L, nrow = nrow(L), ncol = ncol(L))
  for (i in 1:nrow(L)){
    for (j in 1:nrow(L)){
      if (i != j) R[i, j] = abs(sum((Z[,i] - Z[,j])^2))
    }
  }
  return(R)
}


#' @name EffRPar
#' @title Effective Resistance Approximation with RAM Features and Parallel Processing
#' @description  Approximately calculate the effective resistances between all nodes of a graph. \cr
#' Based on algorithm from Spielman-Srivastava (2011) \cr
#' This version of the effective resistance calculator is more RAM efficient and uses the parallel processing package **snow**
#' @import igraph cPCG
# Input:
#' @param network
#' An adjacency matrix or edge list of the network; edge list should be in the format |in|out|weight|
#' @param epsilon
#' Parameter to control resolution of the approximation
#' @param workers
#' The number of parallel processes to run
#' @param n
#' The number of nodes; default is \code{NULL} and should not be changed unless the network is disconnected
# Output:
#' @return A matrix, \code{R}, such that \code{R_ij} gives the approximate effective resistance between node \code{i} and \code{j}
#' @examples
#' A <- ER_gen(n=100, p=0.1, weights=1)
#' R <- EffRPar(network=A, epsilon=0.1, workers=16)
#' @export
EffRPar <- function(network, epsilon, workers=strtoi(Sys.getenv('NUMBER_OF_PROCESSORS')), n = NULL){
  if (nrow(network) != ncol(network)){
    if (is.null(n)) n = max(network[,1:2])
    network = EList_Mtrx(network, n)
  } else {
    n = nrow(network)
  }
  m = dim(Mtrx_EList(network))[1]
  tolProb = 0.5
  L = as.matrix(diag(rowSums(abs(network))) - network)
  B = sVIM(network) #Could be faster...
  W = WDiag(network)
  k = ceiling((log(n,2)) / (epsilon))
  batch = ceiling(k/workers)

  #batcher = function(x){
  #  t(as.matrix(x%*%as.matrix(W)%*%as.matrix(B)))
  #}

  solver = function(x){
    cPCG::cgsolve(L, x)
  }

  cc = snow::makeSOCKcluster(workers)
  snow::clusterExport(cc, 'L', envir=environment())
  snow::clusterExport(cc, 'cgsolve', envir=environment())#Needed?

  Z = matrix(NA, nrow=1, ncol=batch*workers*n)
  x = 1
  y = workers * n
  for (i in 1:batch){
    Ys = matrix(0L, nrow=workers, ncol=n)
    for (j in 1:workers){
      ons = stats::runif(m, 0, 1) > tolProb
      ons = ons - !(ons)
      ons = ons/sqrt(k)
      Y_r = t(as.matrix(ons%*%W%*%B))
      Ys[j, ] = Y_r
    }
    #ons = runif(m*workers, 0, 1) > tolProb
    #ons = ons - !(ons)
    #ons = ons/sqrt(k)
    #ons = matrix(ons, nrow=m, ncol=workers)
    #Ys = snow::parRapply(cc, ons, batcher)

    Zc_s = snow::parRapply(cc, Ys, solver)

    Z[,x:y] = Zc_s #May take a long time | could use apply
    x = x + (workers * n)
    y = y + (workers * n)
  }
  snow::stopCluster(cc)
  rm(cc)
  Z = matrix(Z, nrow=workers*batch, ncol=n, byrow=TRUE)
  R = matrix(0L, nrow = nrow(L), ncol = ncol(L))
  for (i in 1:nrow(L)){
    for (j in 1:nrow(L)){
      if (i != j) R[i, j] = abs(sum((Z[,i] - Z[,j])^2))
    }
  }

  return(R)
}


#' @name normprobs
#' @title Normalize Probabilities
#' @description  Normalize a list of probabilities such that \code{sum(probs)=1}
# Input:
#' @param P
#' List of probabilities
# Output:
#' @return A list of probabilities, \code{P_n}, such that they sum to 1
#' @export
normprobs <- function(P){
  p_f = 1/ sum(P)
  P_n = c()
  for (p in P){
    P_n = append(P_n, p_f * p)
  }

  return(P_n)
}

#' @name EffR_List
#' @title Effective Resistance Matrix to Effective Resistance Edge List
#' @description  Create a vector of effective resistance edge list from the output of \link[simplifyNet]{EffR}, or other versions
# Input:
#' @param R
#' Matrix of effective resistance values
#' @param network
#' Adjacency matrix of the corresponding network
#' @param n
#' The number of nodes; default is \code{NULL} and should not be changed unless the network is disconnected
# Output:
#' @return Vector, \code{R_v}, of effective resistance values
#' @export
EffR_List <- function(R, network, n = NULL){
  #if (nrow(network) != ncol(network)){
  #  A = EList_Mtrx(network, n = n)
  #}
  R_v = c()
  network[upper.tri(network)] = 0
  R[upper.tri(R)] = 0
  for (i in 1:nrow(network)){
    for (j in 1:nrow(network)){
      if (network[i,j] > 0){
        R_v = append(R_v, R[i,j])
      }
    }
  }

  return(R_v)
}


#' @name EffRSparse
#' @title Effective Resistance Spectral Sparsification
#' @description Create a effective resistance sparsifer \cr
#' From Spielman and Srivastava 2008 \cr
#' Effective resistance edge values are used to randomly sample edges with replacement. After \code{O(nlog^c(n)/epsilon^2)} samples with replacement, the algorithm produces a weighted
#' subgraph such that the Laplacian is approximately multiplicity preserved within +/-\code{epsilon}.
# Input:
#' @param network
#' Adjacency matrix of the given network
#' @param q
#' Number of samples to independently sample with replacement
#' @param R
#' Matrix of effective resistances; see \link[simplifyNet]{EffR}
#' @param n
#' The number of nodes; default is \code{NULL} and should not be changed unless the network is disconnected
# Output:
#' @return Generate an effective resistance spectral sparsifer adjacency matrix, \code{H}
#' @examples
#' A <- ER_gen(n=100, p=0.1, weights=1)
#' R <- EffRPar(network=A, epsilon=0.1, workers=16)
#' m = length(Mtrx_EList(A))
#' q = ceiling((m*(log(m))^1)/(0.1^2))
#' G <- EffRSparse(network=A, q=q, R=R)
#' @export
EffRSparse <- function(network, q, R, n = NULL){
  if (nrow(network) == ncol(network)){
    n = nrow(network)
    A = Mtrx_EList(network)
  } else {
    A = network
  }
  W_list = A[,3]
  R_list = EffR_List(R, EList_Mtrx(A, n))
  P = c()
  for (i in 1:nrow(A)){
    w_e = W_list[i]
    R_e = R_list[i]
    P = append(P, (w_e * R_e) / (n - 1))
  }
  P_n = normprobs(P)
  A = cbind(A, P_n)
  C = A[sample(1:nrow(A), size = q, replace = TRUE, prob = P_n),]
  if(is.null(nrow(C))){
    C = matrix(C, nrow = 1, ncol = 4, byrow=TRUE)
  }
  H = matrix(0L, nrow = n, ncol = n)
  for (j in 1:q){
    w_e = C[j, 3]
    p_e = C[j, 4]
    H[C[j,1], C[j,2]] = (w_e / (q * p_e)) + H[C[j,1], C[j,2]]
  }
  H = H + t(H)

  return(H)
}

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
#' @param  m
#' Number edges in the given network
#' @param epsilon
#' Used to specify resolution of the approximation with regard to the Johnson-Lindenstrauss lemma
# Output:
#' @return Generates a random +/- 1/sqrt(k) Bernoulli matrix, \code{Q}
#' @export
QGen <- function(m, epsilon){
  k = ((24 * log(m))) / (epsilon ^ 2)
  k_d = ceiling(k)
  x = 1 / sqrt(k)
  Q = matrix(sample(c(-x,x), size = k_d * m, replace = TRUE),
             nrow = k_d,
             ncol = m)

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
#' @export
EffR <- function(network, epsilon, n = NULL){
  if (nrow(network) != ncol(network)){
    if (is.null(n)) n = max(network[,1:2])
    network = EList_Mtrx(network, n)
  }
  A = methods::as(network, 'matrix')
  B = sVIM(network)
  L = methods::as(igraph::laplacian_matrix(graph_from_adjacency_matrix(A, mode = 'undirected', weighted = T), normalized = FALSE), 'matrix')
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
  R = matrix(0L, nrow = nrow(A), ncol = ncol(A))
  for (i in 1:nrow(A)){
    for (j in 1:nrow(A)){
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
#' @export
EffR2 <- function(network, epsilon, n = NULL){
  if (nrow(network) != ncol(network)){
    if (is.null(n)) n = max(network[,1:2])
    network = EList_Mtrx(network, n)
  }
  A = methods::as(network, 'matrix')
  B = sVIM(network)
  L = methods::as(igraph::laplacian_matrix(graph_from_adjacency_matrix(A, mode = 'undirected', weighted = T), normalized = FALSE), 'matrix')
  W = WDiag(network)
  WB = sqrt(W) %*% B
  edge_num = nrow(B)
  node_num = nrow(A)
  k = ((24 * log(edge_num))) / (epsilon ^ 2)
  k_D = ceiling(k)
  x = 1 / sqrt(k)
  rm(W, B)
  Z = matrix(NA, nrow = k_D, ncol = node_num)
  for (i in 1:k_D){
    Y_r = matrix(0L, nrow=1, ncol=node_num)
    C = matrix(sample(c(-x,x), size = edge_num, replace = TRUE), nrow=1)
    for (j in 1:node_num){
      WB_c = matrix(WB[ ,j])
      entry = C %*% WB_c
      Y_r[1,j] = entry
    }
    Z[i,] = cPCG::cgsolve(L, Y_r)
  }
  R = matrix(0L, nrow = nrow(A), ncol = ncol(A))
  for (i in 1:nrow(A)){
    for (j in 1:nrow(A)){
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
#' @export
EffRPar <- function(network, epsilon, workers=strtoi(Sys.getenv('NUMBER_OF_PROCESSORS')), n = NULL){
  if (nrow(network) != ncol(network)){
    if (is.null(n)) n = max(network[,1:2])
    network = EList_Mtrx(network, n)
  }
  A = methods::as(network, 'matrix')
  B = sVIM(network)
  L = methods::as(igraph::laplacian_matrix(graph_from_adjacency_matrix(A, mode = 'undirected', weighted = T), normalized = FALSE), 'matrix')
  W = WDiag(network)
  WB = sqrt(W) %*% B
  edge_num = nrow(B)
  node_num = nrow(A)
  k = ((24 * log(edge_num))) / (epsilon ^ 2)
  k_D = ceiling(k)
  batch = ceiling(k_D/workers)
  x = 1 / sqrt(k)
  rm(W, B)
  Z = matrix(NA, nrow = k_D, ncol = node_num)

  solver = function(x){
    cPCG::cgsolve(L, x)
  }

  cc = snow::makeSOCKcluster(workers)
  snow::clusterExport(cc, 'L', envir=environment())
  snow::clusterExport(cc, 'cgsolve', envir=environment()) #Needed?

  Z = c()#Pre allocate vector
  for (i in 1:batch){
    Ys.df = data.frame(matrix(NA, nrow=workers, ncol=node_num))
    for (j in 1:workers){
      Y_r = matrix(0L, nrow=1, ncol=node_num)
      C = matrix(sample(c(-x,x), size = edge_num, replace = TRUE), nrow=1)
      for (k in 1:node_num){
        WB_c = matrix(WB[ ,k])
        entry = C %*% WB_c
        Y_r[1,k] = entry #apply statement
      }
      Ys.df[j, ] = Y_r
    }

    Zc_s = matrix(snow::parRapply(cc, Ys.df, solver), nrow=workers, ncol=node_num)

    Z = append(Z, Zc_s) #May take a long time | could use apply
  }
  snow::stopCluster(cc)
  Z = matrix(Z, nrow=workers*batch, ncol=node_num, byrow=TRUE)
  R = matrix(0L, nrow = nrow(A), ncol = ncol(A))
  for (i in 1:nrow(A)){
    for (j in 1:nrow(A)){
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
  A[upper.tri(A)] = 0
  R[upper.tri(R)] = 0
  for (i in 1:nrow(A)){
    for (j in 1:nrow(A)){
      if (A[i,j] > 0){
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
#' @export
EffRSparse <- function(network, q, R, n = NULL){
  if (nrow(network) == ncol(network)){
    n = nrow(network)
    A = Mtrx_EList(network)
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
  H = matrix(0L, nrow = n, ncol = n)
  for (j in 1:q){
    w_e = C[j, 3]
    p_e = C[j, 4]
    H[C[j,1], C[j,2]] = (w_e / (q * p_e)) + H[C[j,1], C[j,2]]
  }
  H = H + t(H)

  return(H)
}

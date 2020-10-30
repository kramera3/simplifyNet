# Author: Alexander M. Mercier

# From Spielman and Srivastava 2011
library(igraph)
library(cPCG)

source('Matrix_Mod.R')


# Compute random \pm 1\sqrt(k) Bernoulli matrix where k=24log(n)/epsilon^2
# Used to fulfill Johnson-Lindenstrauss lemma
# From Spielman and Srivastava 2011
# Input:
# m - Number of edges
# epsilon - List of weights (default none)
# Output:
# Q - Random \pm 1/sqrt(k) Bernoulli matrix
QGen <- function(m, epsilon){
  k = ((24 * log(m))) / (epsilon ^ 2)
  k_d = ceiling(k)
  x = 1 / sqrt(k)
  Q = matrix(sample(c(-x,x), size = k_d * m, replace = TRUE),
             nrow = k_d,
             ncol = m)
  #Q = as(Q, 'sparseMatrix') #Sparse is larger than non sparse matrix

  return(Q)
}


# Approximately find effective resistances for edges of a graph
# Based on algorithm from Spielman-Srivastava (2011)
# Input:
# network - adj matrix or edge list
# epsilon - parameter to control accuracy
# Output:
# T - matrix of effective resistances for the edges in elist
EffRPar <- function(network, epsilon, workers=strtoi(Sys.getenv('NUMBER_OF_PROCESSORS')), n = NULL){
  if (nrow(network) != ncol(network)){
    if (is.null(n)) n = max(network[,1:2])
    network = Elist_Mtrx(network, n)
  }
  A = as(network, 'matrix')
  B = sVIM(network)
  L = as(igraph::laplacian_matrix(graph_from_adjacency_matrix(A, mode = 'undirected', weighted = T), normalized = FALSE), 'matrix')
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

  Z = c()
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

    Z = append(Z, Zc_s)
  }
  stopCluster(cc)
  Z = matrix(Z, nrow=workers*batch, ncol=node_num, byrow=TRUE)
  R = matrix(0L, nrow = nrow(A), ncol = ncol(A))
  for (i in 1:nrow(A)){
    for (j in 1:nrow(A)){
      if (i != j) R[i, j] = abs(sum((Z[,i] - Z[,j])^2))
    }
  }

  return(R)
}



# Create a effective resistance sparsifer
# From Spielman and Srivastava 2008
# Input:
# adj - Adj matrix
# q - number of samples
# R - Matrix of effective resistances #(other types of edge importance in the future, possibly)#
# Output:
# H - effective resistance sparsifer adj matrix
EffRSparse <- function(network, q, R, n = NULL){
  if (nrow(network) == ncol(network)){
    n = nrow(network)
    A = Mtrx_Elist(network)
  }
  W_list = A[,3]
  R_list = EffR_List(R, Elist_Mtrx(A, n))
  P = c()
  for (i in 1:nrow(A)){
    w_e = W_list[i]
    R_e = R_list[i]
    P = append(P, (w_e * R_e) / (n - 1))
  }
  P_n = normprobs(P)
  A = cbind(A, P_n)
  C = A[sample(nrow(A), size = q, replace = TRUE, prob = P_n),]
  H = matrix(0L, nrow = n, ncol = n)
  for (j in 1:q){
    w_e = C[j, 3]
    p_e = C[j, 4]
    H[C[j,1], C[j,2]] = (w_e / (q * p_e)) + H[C[j,1], C[j,2]]
  }
  H = H + t(H)

  return(H)
}


















}

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

  return(Q)
}


# Approximately find effective resistances for edges of a graph
# Based on algorithm from Spielman-Srivastava (2011)
# Input:
# network - adj matrix or edge list
# epsilon - parameter to control accuracy
# Output:
# T - matrix of effective resistances for the edges in elist
EffR <- function(network, epsilon, n = NULL){
  if (nrow(network) != ncol(network)){
    if (is.null(n)) n = max(network[,1:2])
    network = Elist_Mtrx(network, n)
  }
  A = as(network, 'sparseMatrix')
  B = sVIM(network)
  L = as(laplacian_matrix(graph_from_adjacency_matrix(A), normalized = FALSE), 'matrix')
  W = WDiag(network)
  Q = QGen(nrow(B), epsilon)
  rm(Q, W, B)
  Y = as((Q %*% (sqrt(W)) %*% B), 'matrix')
  Z = matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
  for (i in 1:nrow(Y)){
    Y_v = matrix(Y[i, ])
    Z[i,] = cgsolve(L, Y_v)
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


# Normalize probs such that sum(probs)=1
# Input:
# P - list of probs
# Output:
# P_n - list of probs that sum to 1
normprobs <- function(P){
  p_f = 1/ sum(P)
  P_n = c()
  for (p in P){
    P_n = append(P_n, p_f * p)
  }

  return(P_n)
}


# Create a list of edge R_eff
# Input:
# R - matrix of R_effs
# A - adj matrix
# Output:
# R_v - vector of edge R_eff
EffR_List <- function(R, network, n = NULL){
  #if (nrow(network) != ncol(network)){
  #  A = Elist_Mtrx(network, n = n)
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

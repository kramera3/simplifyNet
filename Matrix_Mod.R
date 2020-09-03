# Author: Alexander M. Mercier

# Misc. network matrix modifications
library(Matrix)


# Adj matrix to edge list
# Input:
# A - Adj matrix
# Output:
# E_List - Matrix of edges and edge weights
Mtrx_Elist <- function(A){
  n = dim(A)
  A[lower.tri(A, diag=TRUE)] = 0
  edges = which(A != 0, arr.ind = TRUE)
  weights = A[ which(!A == 0)]
  E_List = cbind(edges, weights)
  colnames(E_List) <- c('row', 'col', 'weight')

  return(E_List)
}


# Edge list to adj matrix
# Input:
# E_List - Adj matrix where edge weights are one
# n - number of nodes
# Output:
# E_List - List of lists giving edges
Elist_Mtrx <- function(E_List, n){
  A = matrix(0L, nrow = n, ncol = n)
  for (i in 1:dim(E_List)[1]){
    n1 = E_List[i, 1]
    n2 = E_List[i, 2]
    A[n1, n2] = E_List[i, 3]
  }
  A = A + t(A)

  return(A)
}

# Compute weight diagonal matrix
# Input:
# network - adj matrix or edge list with weights
# Output:
# W - matrix with weights on the diagonal
WDiag <- function(network){
  if (dim(network)[1] == dim(network)[2]){
    weights = Mtrx_Elist(network)[,'weight']
  } else {
    weights = network[,'weight']
  }
  W = sparseMatrix(i = c(1:length(weights)), j = c(1:length(weights)), x = c(weights))

  return(W)
}


# Compute vertex incidence matrix
# Input:
# network - adj matrix or edge list with weights
# Output:
# B - vertex incidence matrix made up of vertex incidence vectors
sVIM <- function(network){
  if (dim(network)[1] == dim(network)[2]){
    E_List = Mtrx_Elist(network)
  } else {
    E_List = network
  }
  data = c()
  row = c()
  col = c()
  for (i in 1:dim(E_List)[1]){
    p = runif(1, min = 0, max = 1)
    edge = E_List[i, 1:2]
    col = append(col, c(edge[1], edge[2]))
    row = append(row, c(i, i))
    if (p < 0.5){ # Randomly determine head/tail orientation
      data = append(data, c(1,-1))
    } else{
      data = append(data, c(-1, 1))
    }
  }
  B = sparseMatrix(i = row, j = col, x = data)

  return(B)
}






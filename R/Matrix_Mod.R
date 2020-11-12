#' @author  Author: Alexander M. Mercier
#' @rdname Matrix_Mod

# Misc. network matrix modifications

#' @name Mtrx_EList
#' @title Adjacency matrix to edge list
#' @description  Convert an adjacency matrix to an edge list
# Input:
#' @param A
#' Adjacency matrix
# Output:
#' @return An edge list \code{E_List} of adjacency matrix \code{A}.
#' @export
Mtrx_EList <- function(A){
  n = dim(A)
  A[lower.tri(A, diag=TRUE)] = 0
  edges = which(A != 0, arr.ind = TRUE)
  weights = A[which(!A == 0)]
  E_List = cbind(edges, weights)
  colnames(E_List) <- c('row', 'col', 'weight')

  return(E_List)
}

#' @name EList_Mtrx
#' @title Edge list to adjacency matrix
#' @description  Convert an edge list to an adjacency matrix
# Input:
#' @param E_List
#' Edge list formatted |in|out|weight|
#' @param n
#' Number of nodes
# Output:
#' @return Adjacency matrix constructed from edge list \code{E_List}
#' @export
EList_Mtrx <- function(E_List, n){
  A = matrix(0L, nrow = n, ncol = n)
  for (i in 1:dim(E_List)[1]){
    n1 = E_List[i, 1]
    n2 = E_List[i, 2]
    A[n1, n2] = E_List[i, 3]
  }
  A = A + t(A)

  return(A)
}


#' @name WDiag
#' @title Weight Diagonal Matrix
#' @description  Compute the weight diagonal matrix of a network
# Input:
#' @param network
#' Adjacency matrix or edge list of a network; edge list should be formatted |in|out|weight|
# Output:
#' @return Matrix, \code{W}, with network weights on the diagonal
#' @export
WDiag <- function(network){
  if (dim(network)[1] == dim(network)[2]){
    weights = Mtrx_EList(network)[,'weight']
  } else {
    weights = network[,'weight']
  }
  W = Matrix::sparseMatrix(i = c(1:length(weights)), j = c(1:length(weights)), x = c(weights))

  return(W)
}


#' @name sVIM
#' @title Vertex Incidence Matrix
#' @description  Compute vertex incidence matrix of a given network
# Input:
#' @param network
#' Adjacency matrix or edge list with weights
# Output:
#' @return A vertex incidence matrix, \code{B},  made up of vertex incidence vectors
#' @export
sVIM <- function(network){
  if (dim(network)[1] == dim(network)[2]){
    E_List = Mtrx_EList(network)
  } else {
    E_List = network
  }
  data = c()
  row = c()
  col = c()
  for (i in 1:dim(E_List)[1]){
    p = stats::runif(1, min = 0, max = 1)
    edge = E_List[i, 1:2]
    col = append(col, c(edge[1], edge[2]))
    row = append(row, c(i, i))
    if (p < 0.5){ # Randomly determine head/tail orientation
      data = append(data, c(1,-1))
    } else{
      data = append(data, c(-1, 1))
    }
  }
  B = Matrix::sparseMatrix(i = row, j = col, x = data)

  return(B)
}






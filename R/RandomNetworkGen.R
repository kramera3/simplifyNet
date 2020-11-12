#' @author Author: Alexander M. Mercier
#' @rdname RandomNetworkGen

# Create a variety of random networks, saving as either edge list or adj matrix

#' @name ER_gen
#' @title Erdos-Renyi Random Network Generation
#' @description Generate an Erdos-Renyi model with specified edge weights and edge weight probabilities
# Input:
#' @param n
#' The number of nodes in the Erdos-Renyi model
#' @param p
#' Probability of linking two nodes with an edge
#' @param weights
#' A vector of edge weights for Erdos-Renyi model to sample from
#' @param w.prob
#' Vector of probabilities to sample edge weights with
#' @param self.edge
#' To include self edges within the Erdos-Renyi model
# Output:
#' @return An adjacency matrix, \code{A}, of the Erdos-Renyi model
#' @export
ER_gen <- function(n, p, weights, w.prob = NULL, self.edge = FALSE){
  #set.seed(001)
  A = matrix(0L, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      if (i <= j) {
        theta = stats::runif(1, min = 0, max = 1)
        if (theta < p){
          if (self.edge == TRUE){
            if (i == j) {A[i, j] = 2 * sample(x = weights, size = 1, replace = TRUE, prob = w.prob)}
            else {A[i, j] = sample(x = weights, size = 1, replace = TRUE, prob = w.prob)}
          } else if (i != j){
            A[i, j] = sample(x = weights, size = 1, replace = TRUE, prob = w.prob)
          }
        }
      }
    }
  }
  A = A + t(A)

  return(A)
}


#' @name ConfigNetwork_gen
#' @title Configuration Model Generation
#' @description  Generate a configuration model from a given degree sequence
# Input:
#' @param k
#' Degree list \cr
#' ##The degrees must sum to an even number
#' @param weights
#' A vector of edge weights for configuration model to sample from
#' @param w.prob
#' A vector of edge weight probabilities to sample edge weights with
# Output:
#' @return An adjacency matrix, \code{A}, of configuration network
#' @export
ConfigNetwork_gen <- function(k, weights, w.prob = NULL){
  if (Reduce('+', k) %% 2 != 0){stop("Please use a degree distribution that sums to an even number")} # The sum of the degrees must be even
  L = list()
  A = matrix(0L, nrow = length(k), ncol = length(k))
  for (i in 1:length(k)){
    for (j in 1:k[[i]]){
      L = append(L, i) # Create stub list
    }
  }
  set.seed
  L = sample(L) # Randomly shuffle stublist
  while (length(L) > 0) {
    n1 = L[[1]]
    L[1] = NULL
    n2 = L[[1]]
    L[1] = NULL
    if (n1 == n2){
      A[n1, n2] = (2 * sample(x = weights, size = 1, replace = TRUE, prob = w.prob)) + A[n1, n2] # Self edges have an added degree of two
    } else {
      A[n1, n2] = A[n2, n1] = sample(x = weights, size = 1, replace = TRUE, prob = w.prob) + A[n1, n2]
    }
  }

  return(A)
}


# Create a random partition of nodes with given community sizes
# Input:
# L - node list
# s - list of community sizes
# Output:
# P - list of lists; a partition based on L and s
l_g_partition <- function(L, s){
  P = list()
  cs = 1
  for (i in 1:length(s)){
    C = list()
    for (j in 1:s[[i]]){
      n = cs
      #L[[1]][j] = NULL
      C = append(C, n) # Add community partition to a list of partitions
      cs = cs + 1
    }
    P = append(P, list(C))
  }

  return(P)
}


#' @name SBM_gen
#' @title Stochastic Block Model Generation
#' @description  Generate a stochastic block model from specified community sizes, between community probability of an edge, and within community probability of an edge
# Input:
#' @param s
#' List of community sizes
#' @param W
#' Matrix of probabilities such that \code{W_ij} is the probability of connecting community \code{i} and \code{j} with an edge
#' @param weights
#' Vector of edge weights for the stochastic block model to sample from
#' @param w.prob
#' Vector of probabilities to sample edge weights with
#' @param self.edge
#' To include self edges within the Erdos-Renyi model
# Output:
#' @return An adjacency matrix, \code{A}, of a stochastic block model
#' @export
SBM_gen <- function(s, W, weights, w.prob = NULL, self.edge = FALSE){
  n = Reduce('+', s) # Find total number of nodes
  A = matrix(0L, nrow = n, ncol = n)
  nodelist = list(1:n) # Create a node list
  P = l_g_partition(nodelist, s) # Create partition of node list
  for (i in 1:length(P)){
    for (j in 1:length(P)){
      C1 = P[[i]]
      C2 = P[[j]]
      for (x in 1:length(C1)){
        for (y in 1:length(C2)){
          if (x <= y){
            p = W[i, j]
            theta = stats::runif(1, min = 0, max = 1)
            if (theta < p){
              n = C1[[x]]
              m = C2[[y]]
              if (n != m){
                A[n, m] = sample(x = weights, size = 1, replace = TRUE, prob = w.prob) + A[n, m]
                if (i == j){
                  A[n, m] = 0 + A[n, m]
                }
              } else if (self.edge == TRUE) {
                A[n, n] = (2 * sample(x = weights, size = 1, replace = TRUE, prob = w.prob)) + A[n, n]
                if (i == j){
                  A[n, m] = 0 + A[n, m]
                }
              }
            }
          }
        }
      }
    }
  }
  A = A + t(A)

  return(A)
}


#' @name RanGeo_gen
#' @title Random Geometric Network Generation
#' @description  Generate a random geometric network with specified maximum edge distance
#' @import fields
# Input:
#' @param n
#' The number of nodes in the random geometric network
#' @param r
#' Maximum euclidean distance for adding an edge in the metric space \code{[0,1)}
#' @param weights
#' Edge weights for random geometric network to sample from
#' @param w.prob
#' Vector of probabilities to sample edge weights with
# Output:
#' @return The adjacency matrix, \code{A}, of the random matrix with a data.frame, \code{pos}, of the x,y coordinates for each node
#' @export
RanGeo_gen <- function(n, r, weights, w.prob = NULL){
  A = matrix(0L, nrow = n, ncol = n)
  pos = data.frame(matrix(NA, nrow = n, ncol = 2))
  pos$X1 = stats::runif(n, min = 0, max = 1)
  pos$X2 = stats::runif(n, min = 0, max = 1)
  dist_m = fields::rdist(x1 = pos[,1:2])
  for (i in 1:n){
    for (j in 1:n){
      if (i < j) {
        if (dist_m[i,j] <= r){
          A[i, j] = A[j, i] = sample(x = weights, size = 1, replace = TRUE, prob = w.prob)
        }
      }
    }
  }

  return(list(A,pos))
}



# Author: Alexander M. Mercier

# Create a variety of random networks, saving as either edge list or adj matrix

library(tidyr)
library(fields)

# Generate an Erdos-Renyi Model
# Input:
# n - the number of nodes
# p - probability of adding an edge
# weights - edge weights for ER network to sample from
# w.prob - list of probabilities to sample edge weights with
# Output:
# A or E - adj matrix or edge list of ER network
ER_gen <- function(n, p, weights, w.prob = NULL, self_edge = FALSE){
  #set.seed(001)
  A = matrix(0L, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      if (i <= j) {
        theta = runif(1, min = 0, max = 1)
        if (theta < p){
          if (self_edge == TRUE){
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


# Generate a Configuration Model
# Input:
# k - degree list
# weights - edge weights for ER network to sample from
# w.prob - list of probabilities to sample edge weights with
# Output:
# A - adj matrix of Configuration network
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


# Create a SBM
# Input:
# s - list of community sizes
# W - matrix of probabilities where W_ij is the prob of connecting community i and j
# weights - edge weights for ER network to sample from
# w.prob - list of probabilities to sample edge weights with
# Output:
# A - adj matrix of SBM network
SBM_gen <- function(s, W, weights, w.prob = NULL, self_edge = FALSE){
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
            theta = runif(1, min = 0, max = 1)
            if (theta < p){
              n = C1[[x]]
              m = C2[[y]]
              if (n != m){
                A[n, m] = sample(x = weights, size = 1, replace = TRUE, prob = w.prob) + A[n, m]
                if (i == j){
                  A[n, m] = 0 + A[n, m]
                }
              } else if (self_edge == TRUE) {
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


# Generate random undirected geometric graph
# Input:
# n - the number of nodes
# r - minimum dist to add edge
# weights - edge weights for ER network to sample from
# w.prob - list of probabilities to sample edge weights with
# Output:
# A - the adj matrix of the graph
# pos - data.frame of x, y coordinate for each node
RanGeo_gen <- function(n, r, weights, w.prob = NULL){
  A = matrix(0L, nrow = n, ncol = n)
  pos = data.frame(matrix(NA, nrow = n, ncol = 2))
  pos$X1 = runif(n, min = 0, max = 1)
  pos$X2 = runif(n, min = 0, max = 1)
  dist_m = rdist(x1 = pos[,1:2])
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




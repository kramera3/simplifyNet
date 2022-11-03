#' @name net.as
#' @title Network format converter
#' @description  Convert a network in weighted adjacency matrix, edge list, or igraph to a weighted adjacency matrix, edge list, or igraph format.
# Input:
#' @param network
#' Weighted adjacency matrix, weighted \code{igraph} network, or edge list formatted | n1 | n2 | weight | with colnames \code{c("n1", "n2", "weight")}.
#' @param net.to
#' Specifics to what format the imputed network is to be converted: \cr
#' (1) 'E_List' convert to an edge list of the format | n1 | n2 | weight | with colnames \code{c("n1", "n2", "weight")}.\cr
#' (2) 'Adj' convert to a weighted adjacency matrix.\cr
#' (3) 'igraph' convert to a weighted igraph object.
#' @param directed
#' If \code{TRUE}, specifies that the inputted network is directed. Default is \code{FALSE}.
# Output:
#' @return A network of the format specified by \code{net.to}.
#' @author Alexander Mercier
#' @export
net.as <- function(network, net.to = "E_List", directed = FALSE){
  if(net.to == "E_List"){
    if(any(class(network) %in% c("matrix", "Matrix", "data.frame")) &&
       dim(network)[1] == dim(network)[2] &&
       !(identical(colnames(network), c("n1", "n2", "weight")))){
      net = simplifyNet::Mtrx_EList(network, directed = directed)
    } else if(any(class(network) == "igraph")){
      edges = apply(igraph::as_edgelist(network), 2, as.numeric)
      net = cbind(edges, igraph::E(network)$weight)
      colnames(net) <- c("n1", "n2", "weight")
    } else {net = network}
  } else if(net.to == "Adj"){
    if(any(class(network) %in% c("matrix", "Matrix", "data.frame")) &&
       identical(colnames(network), c("n1", "n2", "weight"))){
      net = simplifyNet::EList_Mtrx(network, directed = directed)
    } else if(any(class(network) == "igraph")){
      net = igraph::as_adjacency_matrix(network, attr="weight", sparse=TRUE)
    } else {net = network}
  } else if(net.to == "igraph"){
    if(any(class(network) %in% c("matrix", "Matrix", "data.frame")) &&
       dim(network)[1] == dim(network)[2] &&
       !(identical(colnames(network), c("n1", "n2", "weight")))){
      if(directed){
        mode = "directed"
      } else {mode = "undirected"}
      net = igraph::graph_from_adjacency_matrix(network, mode = mode, weighted = TRUE)
    } else if(any(class(network) %in% c("matrix", "Matrix", "data.frame")) &&
              identical(colnames(network), c("n1", "n2", "weight"))){
      df <- data.frame(network)
      net <- igraph::graph.data.frame(df, directed = directed)
    } else {net = network}
  }
  return(net)
}

#' @name igraph.to.elist
#' @title igraph to EList
#' @description Convert igraph object to edge list.
#' @keywords internal
#' @details Intended as internal function.
#' @return Return an edge list with weights from an igraph object.
#' @author Alexander Mercier
#' @export
igraph.to.elist <- function(g){
  if(is.null(igraph::E(g)$weight)){
    E_List = igraph::as_edgelist(g)
    E_List = cbind(E_List, 1L) #If not weights, add edge weights 1
  } else {
    E_List = cbind(igraph::as_edgelist(g), igraph::E(g)$weight)
  }
  return(E_List)
}

#' @name Mtrx_EList
#' @title Adjacency matrix to edge list
#' @description  Convert an adjacency matrix to an edge list.
# Input:
#' @param A
#' Weighted adjacency matrix.
#' @param directed
#' Specifies if the network is directed or undirected. Default is set to undirected.
# Output:
#' @return An edge list, \code{E_List}, of adjacency matrix, \code{A}, of the form | n1 | n2 | weight |.
#' @author Alexander Mercier
#' @export
Mtrx_EList <- function(A , directed=FALSE){
  n = dim(A)
  if(directed == FALSE){
    A[upper.tri(A, diag=TRUE)] = 0 #Functions below only work due to col, row order
  }
  edges = which(methods::as(A != 0, "matrix"), arr.ind = TRUE) #Better way of doing this?
  weights = A[which(methods::as(!A == 0, 'matrix'))]
  E_List = cbind(edges, weights)
  colnames(E_List) <- c('n1', 'n2', 'weight')

  return(E_List)
}


#' @name EList_Mtrx
#' @title Edge list to adjacency matrix
#' @description  Convert an edge list to an adjacency matrix.
# Input:
#' @param E_List
#' Edge list formatted | n1 | n2 | weight |.
#' @param directed
#' Specifies if the network is directed or undirected. Default is set to undirected.
#' @param n
#' Specify number of nodes. Default is \code{max(E_List[,1:2])}.
# Output:
#' @return Adjacency matrix constructed from edge list, \code{E_List}, of the class dgCMatrix.
#' @author Alexander Mercier
#' @export
EList_Mtrx <- function(E_List, directed=FALSE, n){
  if(directed == FALSE){
    edges1 = c(E_List[,1], E_List[,2])
    edges2 = c(E_List[,2], E_List[,1])
    weights = c(E_List[,3], E_List[,3])
    if(missing(n)){
      A = Matrix::sparseMatrix(i=edges1, j=edges2, x=weights)
    } else {
      A = Matrix::sparseMatrix(i=edges1, j=edges2, x=weights, dims = c(n, n))
    }
  }
  if(directed == TRUE){
    if(missing(n)){
      A = Matrix::sparseMatrix(i=E_List[,1], j=E_List[,2], x=E_List[,3])
    } else {
      A = Matrix::sparseMatrix(i=E_List[,1], j=E_List[,2], x=E_List[,3], dims = c(n, n))
    }
  }
  return(A)
}

#' @name Lap
#' @title Graph Laplacian calculator
#' @description  Find the graph Laplacian from a weighted adjacency matrix.
# Input:
#' @param A
#' Weighted adjacency matrix.
# Output:
#' @return Graph Laplacian, \code{L}, of the class dgCMatrix.
#' @keywords internal
#' @details Intended as internal function.
#' @author Alexander Mercier
#' @examples
#' A = matrix(c(0,1,1,1,0,1,1,1,0), 3 ,3)
#' L = Lap(A)
#' @export
Lap <- function(A){
  W = Matrix::sparseMatrix(i=1:ncol(A), j=1:ncol(A), x=Matrix::colSums(abs(A)))
  L = W - A
  return(L)
}

#' @name sVIM
#' @title Signed vertex incidence matrix
#' @description  Calculate the signed vertex incidence matrix.
# Input:
#' @param E_List
#' Edge list formatted | n1 | n2 | weight |.
# Output:
#' @return Return a signed, m by n, vertex incidence matrix, \code{B}.
#' @keywords internal
#' @details Intended as internal function.
#' @author Alexander Mercier
#' @examples
#' E_List = matrix(c(1,1,2,2,3,3,1,1,1), 3, 3)
#' colnames(E_List) <- c("n1", "n2", "weight")
#' B = sVIM(E_List)
#' @export
# Find signed-edge vertex incidence matrix, B mxn
sVIM <- function(E_List){
  m = dim(E_List)[1]

  data = c(rep(1,m), rep(-1,m))
  i = c(1:m,1:m)
  j = c(E_List[,1], E_List[,2])
  B = Matrix::sparseMatrix(i=i, j=j, x=data)# Sum of any row of a sVIM will be 0
  return(B)
}


#' @name WDiag
#' @title Diagonal, m by m, weight matrix
#' @description  Calculate the m by m diagonal weight matrix, \code{W}, for an imputed graph.
# Input:
#' @param weights
#' List of edges corresponding to edge list.
# Output:
#' @return Return diagonal weight matrix, \code{W}.
#' @keywords internal
#' @details Intended as internal function.
#' @author Alexander Mercier
#' @export
WDiag <- function(weights){
  W = Matrix::Diagonal(x=weights)
  return(W)
}


#' @name sqrt_WDiag
#' @title Diagonal, m by m, weight matrix with \code{sqrt(w_e)} on the diagonal.
#' @description  Calculate the m by m diagonal weight matrix, \code{W}, for an inputted graph with \code{sqrt(w_e)} on the diagonal.
# Input:
#' @param weights
#' List of edges corresponding to edge list.
# Output:
#' @return Return diagonal weight matrix, \code{W}, with \eqn{\sqrt(w_e)} on the diagonal.
#' @keywords internal
#' @details Intended as internal function.
#' @author Alexander Mercier
#' @export
sqrt_WDiag <- function(weights){
  weights_sqrt = sqrt(weights)
  W = Matrix::Diagonal(x=weights_sqrt)
  return(W)
}

#Method from Koutis et al.
#Three types: exact (ext), original Spielman/Srivastava (spl), and Koutis et al.(kts)(recommended)


#' @name EffR
#' @title Effective resistances calculator
#' @description Calculate or approximate the effective resistances of an inputted, undirected graph. There are three methods. \cr
#' (1) 'ext' which exactly calculates the effective resistances (WARNING! Not ideal for large graphs).\cr
#' (2) 'spl' which approximates the effective resistances of the inputted graph using the original Spielman-Srivastava algorithm.\cr
#' (3) 'kts' which approximates the effective resistances of the inputted graph using the implementation by Koutis et al. (ideal for large graphs where memory usage is a concern).\cr
#' The relative fidelity of the approximation methods is governed by the variable epsilon.
# Input:
#' @param network
#' Weighted adjacency matrix, weighted \code{igraph} network, or edge list formatted | n1 | n2 | weight | with colnames \code{c("n1", "n2", "weight")}.
#' @param epsilon
#' Variable epsilon governs the relative fidelity of the approximation methods 'spl' and 'kts'. The smaller the value the greater the fidelity of the approximation. Default value is 0.1.
#' @param type
#' There are three methods. \cr
#' (1) 'ext' which exactly calculates the effective resistances (WARNING! Not ideal for large graphs).\cr
#' (2) 'spl' which approximates the effective resistances of the inputted graph using the original Spielman-Srivastava algorithm.\cr
#' (3) 'kts' which approximates the effective resistances of the inputted graph using the implementation by Koutis et al. (ideal for large graphs where memory usage is a concern).
#' @param tol
#' Tolerance for the linear algebra (conjugate gradient) solver to find the effective resistances. Default value is 1e-10.
# Output:
#' @return Return either exact or approximate effective resistances for each edge in the same order as "weight" in the edge list.
#' @author Alexander Mercier
#' @references
#' Spielman, D. A., & Srivastava, N. (2011). Graph sparsification by effective resistances. SIAM Journal on Computing, 40(6), 1913-1926.\cr
#' Koutis, I., Miller, G. L., & Peng, R. (2014). Approaching optimality for solving SDD linear systems. SIAM Journal on Computing, 43(1), 337-354.
#' @details The fidelity of the effective resistance approximation decreases with a decrease in \code{epsilon}. However, it is more important for sparsification by effective resistances that the approximations be roughly equivalent relative to one another, as they will be normalized in a probability distribution where exact values are not needed. \cr
#' The number of edges contained in the sparsifier will be less than or equal to the number of samples taken, \code{q}.
#' @examples
#' E_List = matrix(c(1,1,2,2,3,3,1,1,1), 3, 3) #Triangle graph, \eqn{K_3}, with edge weights equal to 1
#' effR = simplifyNet::EffR(E_List, epsilon = 0.1, type = 'kts', tol = 1e-10)
#' @export
EffR <- function(network, epsilon = 0.1, type = 'kts', tol=1e-10){
  E_List = net.as(network, "E_List")#Always directed
  m = dim(E_List)[1]
  n = max(E_List[,1:2]) #1 indexing

  A = simplifyNet::EList_Mtrx(E_List)
  L = simplifyNet::Lap(A)
  B = simplifyNet::sVIM(E_List)
  W = simplifyNet::sqrt_WDiag(E_List[,3])


  #Exact effR values:

  if(type=='ext'){ #Exact method (ext)
    effR = matrix(0L, 1, m)

    for(i in 1:m){
      Br = matrix(B[i,], 1, m, byrow=TRUE) # row vector
      Z = sanic::solve_cg(L, t(Br), tol=tol)
      R_eff = Br %*%Z
      effR[,i] = R_eff
    }
    return(effR)
  }

  #Original Spielman-Srivastava implementation
  if(type=='spl'){
    scale = ceiling(log2(n))/epsilon
    Q1 = Matrix::rsparsematrix(scale, m, 1, rand.x=stats::runif) > 0.5
    Q2 = Matrix::rsparsematrix(scale, m, 1, rand.x=stats::runif) > 0
    Q_not = Q1-Q2
    Q = Q1 + Q_not #create random -1, 1 matrix
    Q = Q / sqrt(scale)

    SYS = Q%*%W%*%B
    Z = matrix(0L, scale, n)

    for(i in 1:scale){
      SYSr = matrix(SYS[i,], 1, n)
      Z[i, ] = sanic::solve_cg(L, t(SYSr), tol=tol)
    }
    effR = colSums((Z[,E_List[,1]] - Z[,E_List[,2]])^2)
    return(matrix(effR, ncol=m))
  }

  #Koutis et al. implementation
  if(type=='kts'){
    scale = ceiling(log2(n))/epsilon
    effR = matrix(0L, 1, m)

    for(i in 1:scale){
      ones1 = Matrix::rsparsematrix(1, m, 1, rand.x=stats::runif) > 0.5
      ones2 = Matrix::rsparsematrix(1, m, 1, rand.x=stats::runif) > 0
      ones_not = ones1 - ones2
      ones = ones1 + ones_not
      q = ones / sqrt(scale)

      b = q%*%W%*%B

      Z = sanic::solve_cg(L, Matrix::t(b), tol=tol)

      effR = effR + abs((Z[E_List[,1]] - Z[E_List[,2]])^2)

    }
    return(effR)
  }
}


#' @name normProbs
#' @title Normalize probabilities
#' @description Normalize numerics to probabilities such that the sum of the vector equals 1.
# Input:
#' @param P
#' Numerics to be normalized to probabilities.
# Output:
#' @return Return probabilities where their sum is 1.
#' @author Alexander Mercier
#' @keywords internal
#' @details Intended as internal function.
#' @export
normProbs <- function(P){
  prob_frac = 1 / sum(P)
  nP = prob_frac * P
  return(nP)
}


#' @name EffRSparse
#' @author Daniel A. Spielman,
#' @title Sparsification through edge sampling via effective resistances
#' @description Sparsify an undirected network by sampling edges proportional to \eqn{w_e * R_e} where \eqn{w_e} is the weight of edge \eqn{e} and \eqn{R_e} is the effective resistance of edge \eqn{e}.\cr
#' Approximately preserves the graph Laplacian, \code{L}, with increasing fidelity as the number of samples taken increases.
# Input:
#' @param network
#' Weighted adjacency matrix, weighted \code{igraph} network, or edge list formatted | n1 | n2 | weight | with colnames \code{c("n1", "n2", "weight")}.
#' @param q
#' The numbers of samples taken. The fidelity to the original network increases as the number of samples increases, but decreases the sparseness.
#' @param effR
#' Effective resistances corresponding to each edge. Should be in the same order as "weight".
#' @param seed
#' Set the seed to reproduce results of random sampling.
#' @param n
#' The number of nodes in the network. Default is the max node index of the edge list.
# Output:
#' @return A sparsified network, \code{H}, edge list where the number of edges is dependent on the number of samples taken, \code{q}.
#' @author Alexander Mercier
#' @references Spielman, D. A., & Srivastava, N. (2011). Graph sparsification by effective resistances. SIAM Journal on Computing, 40(6), 1913-1926.
#' @details The performance of this method is dependent on the size of the network and fidelity of the effective resistance approximation. The network should be "sufficiently large." \cr
#' For more details, see: https://epubs.siam.org/doi/epdf/10.1137/080734029
#' @examples
#' #Generate random ER graph with uniformly random edge weights
#' g = igraph::erdos.renyi.game(100, 0.1)
#' igraph::E(g)$weight <- runif(length(igraph::E(g)))
#' #Approximate effective resistances
#' effR = simplifyNet::EffR(g)
#' #Use effective resistances to create spectral sparsifier by edge sampling
#' S = simplifyNet::EffRSparse(g, q = 200, effR = effR, seed = 150)
#' sg = simplifyNet::net.as(S, net.to="igraph", directed=FALSE)
#' igraph::ecount(sg)/igraph::ecount(g)#fraction of edges in the sparsifier
#' @export
EffRSparse <- function(network, q, effR, seed, n){
  E_List = simplifyNet::net.as(network, "E_List") #Never should be directed
  if(missing(n)){
    n = max(E_List[,1:2])
  }
  if(!missing(seed)){
    set.seed(seed)
  }
  m = dim(E_List)[1]
  P = matrix(0L, 1, m)
  H_List = matrix(0L, m, 3)
  weights = E_List[,3]

  for(i in 1:m){
    w_e = weights[i]
    R_e = effR[i]
    P[i] = (w_e * R_e) / (n-1)
  }
  nP = simplifyNet::normProbs(P)
  Cdf = data.frame(matrix(0L, m, 5))
  colnames(Cdf) <- c('n1', 'n2', 'weight', 'P', 'i')
  Cdf$n1 = E_List[,1]
  Cdf$n2 = E_List[,2]
  Cdf$weights = E_List[,3]
  Cdf$P = t(nP)
  Cdf$i = 1:m

  C = dplyr::sample_n(Cdf, q, replace=TRUE, weight = nP) #Drop zero rows

  for(i in 1:q){
      n1 = C$n1[i]
      n2 = C$n2[i]
      w_e = C$weights[i]
      p_e = C$P[i]
      index = C$i[i]
      H_List[index, 3] = w_e / (q * p_e) + H_List[index, 3]
      H_List[index, 1:2] = c(n1,n2)
  }
  colnames(H_List) <- c('n1', 'n2', 'weight')
  H_List = H_List[apply(H_List[,-1], 1, function(x) !all(x==0)),]
  return(H_List)
}



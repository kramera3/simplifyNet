#' @author Alexander Mercier
#' @rdname ModelFit

#' @name irefit
#' @title Iterative refitting
#' @description  Iterative sparsifcation based
# Input:
#' @param E_List
#' Edge list of the given network in the format of | node 1 | node 2 | weight |.
#' @param func
#' Model function whose input is the network and whose output is a single real value or a list of reevaluated weights in the first index and a real value in the second.\cr
#' A wrapper function may have to be written.
#' @param tol
#' Allowed error around the original output of \code{func} approximated by the sparsified network within which edges are removed. Specifies if method converges.
#' @param rank
#' Ranking of edges. Lower ranked edges are removed first. Must be the same length as \code{nrow(E_List)}.
#' @param connected
#' If TRUE, connectivity of the network is prioritized over scoring by \code{func}.
#' @param per
#' Percentage of edges to add/remove from the sparsifier at each step.
# Output:
#' @return Sparsified network, \code{func(H)}, which still maintains evaluator function, \code{func}, plus/minus \code{tol}.
#' @export
irefit <- function(E_List, func, tol, rank = 'weight', connected = TRUE, per = 0.5){
  if(rank == 'weight'){
    E_List = E_List[order(E_List[,3], decreasing = FALSE), ]
  } else if(rank == 'random'){
    E_List = E_List[sample(nrow(E_List)), ]
  } else {
    E_List = E_List[order(cbind(E_List, methods::as(rank, 'vector'))[,4], decreasing = FALSE), ]
  }
  org_score = func(E_List)
  if(methods::is(org_score, "list")){
    E_List = cbind(E_List, org_score[[1]])
    org_score = org_score[[2]]
    stepsize = ceiling(per*nrow(E_List))
    index = nrow(E_List) - stepsize
    S_List = E_List[1:index,]
    spr_score = func(S_List)
    S_List = cbind(S_List, spr_score[[1]])
    spr_score = spr_score[[2]]
  }
  stepsize = ceiling(per*nrow(E_List))
  index = nrow(E_List) - stepsize
  S_List = E_List[1:index,]
  spr_score = func(S_List)
  S_List = sparse.step(E_List,
                      S_List,
                      stepsize,
                      spr_score,
                      org_score,
                      func,
                      tol,
                      per,
                      connected)
  return(S_List)
}


sparse.step <- function(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected){
  S = simplifyNet::EList_Mtrx(S_List)
  print(cat(stepsize, nrow(S_List), spr_score, is.connected(S)))
  stepsize = ceiling(per * stepsize)
  if(connected){
    if(abs(spr_score - org_score) > tol || !(is.connected(S))){
      index = nrow(S_List) + stepsize
      S_List = E_List[1:index,]
      spr_score = func(S_List)
      if(methods::is(spr_score, "list")){
        S_List = cbind(S_List, spr_score[[1]])
        spr_score = spr_score[[2]]
      }
      sparse.step(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected)
    } else if(stepsize == 1) {
        return(S_List)
    } else if(abs(spr_score - org_score) <= tol) {
      index = nrow(S_List) - stepsize
      S_List = E_List[1:index,]
      spr_score = func(S_List)
      if(methods::is(spr_score, "list")){
        S_List = cbind(S_List, spr_score[[1]])
        spr_score = spr_score[[2]]
      }
      sparse.step(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected)
    }
  } else {
    if(abs(spr_score - org_score) > tol) {
      index = nrow(S_List) + stepsize
      S_List = E_List[1:index,]
      spr_score = func(S_List)
      if(methods::is(spr_score, "list")){
        S_List = cbind(S_List, spr_score[[1]])
        spr_score = spr_score[[2]]
      }
      sparse.step(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected)
    } else if(abs(spr_score - org_score) <= tol) {
      index = nrow(S_List) - stepsize
      S_List = E_List[1:index,]
      spr_score = func
      if(methods::is(spr_score, "list")){
        S_List = cbind(S_List, spr_score[[1]])
        spr_score = spr_score[[2]]
      }
      sparse.step(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected)
    }
  }
}

#' @name find.neighbors
#' @title Find Neighbors
#' @description  Finds the neighbors of a given node, \code{v}.
# Input:
#' @param Adj
#' Adjacency matrix
#' @param v
#' Node to find the neighbors of
# Output:
#' @return Integers designating node indices of the adjacency matrix for the neighbors of \code{v}
#' @export
find.neighbors <- function(Adj, v){
  neighbors = which(Adj[v, ] != 0, arr.ind = TRUE)
  return(neighbors)
}

#' @name DFS
#' @title Depth First search
#' @description  Iterative depth first search.
# Input:
#' @param Adj
#' Adjacency matrix
#' @param v
#' Node to perform DFS from
#' @param discovered
#' A list of discovered nodes from \code{v}. If initilzing the search, should be all FALSE
# Output:
#' @return Logical of length n where TRUE denotes connected to node \code{v}
#' @export
DFS <- function(Adj, v, discovered){
  discovered[v] <- TRUE
  neighbors <- find.neighbors(Adj, v)
  for(w in neighbors){
    if(!discovered[w]){
      discovered = DFS(Adj, w, discovered)
    }
  }
  return(discovered)
}

#' @name is.connected
#' @title Connectivity of Graph
#' @description  Tests the connectivity of a graph by performing a Depth First Search (DFS) from a random node.
# Input:
#' @param Adj
#' Adjacency matrix
# Output:
#' @return Return TRUE if network is connected and FALSE if not connected.If the network is directed, then this function checks if the network is strongly connected.
#' @export
is.connected <- function(Adj){
  discovered = matrix(FALSE, nrow=1, ncol=ncol(Adj))
  v = sample(1:ncol(Adj), 1)
  discovered = DFS(Adj, v, discovered)
  return(all(discovered))
}

remove.edges <- function(per, E_List, S_List){
  m = nrow(E_List)
  s = nrow(S_List)
  frac = s/m
  new_frac = ceiling(frac * (1-per) * nrow(E_List))
  S_List = E_List[1:new_frac,]
  return(S_List)
}

add.edges <- function(per, E_List, S_List) {
  m = nrow(E_List)
  s = nrow(S_List)
  frac = s/m
  new_frac = ceiling((frac * (1-per) + 1) * nrow(S_List))
  S_List = E_List[1:new_frac,]
  return(S_List)
}

#' @name LeverageScore
#' @title Total statistical leverage score
#' @description  Computes the global statistical leverage score, \eqn{w_e\times R_e}, of a given network
# Input:
#' @param E_List
#' Edge list of the given network in the format of | node 1 | node 2 | weight |
# Output:
#' @return Returns the sum of all statistical leverage scores of each edge, \eqn{\sum_{e \in E} w_e \times R_e}
#' @export
LeverageScore <- function(E_List){
  effR = simplifyNet::EffR(E_List)
  StatLev = effR * E_List[,3]
  return(sum(StatLev))
}

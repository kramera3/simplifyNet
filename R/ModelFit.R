#' @name irefit
#' @title Iterative refitting
#' @description  Iterative sparsifcation based refitting.
# Input:
#' @param network
#' Weighted adjacency matrix, weighted \code{igraph} network, or edge list formatted | n1 | n2 | weight | with colnames \code{c("n1", "n2", "weight")}.
#' @param func
#' Model function whose input is the network and whose output is a single real value or a list of reevaluated weights in the first index and a real value in the second.\cr
#' A wrapper function may have to be written.
#' @param tol
#' Allowed error around the original output of \code{func} approximated by the sparsified network within which edges are removed. Specifies if method converges.
#' @param rank
#' Ranking of edges. Lower ranked edges are removed first. Must be the same length as \code{nrow(E_List)}.
#' @param connected
#' If TRUE, connectivity of the network is prioritized over scoring by \code{func}.
#' @param directed
#' If \code{TRUE}, specifies that the inputted network is directed. Default is \code{FALSE}.
#' @param per
#' Percentage of edges to add/remove from the sparsifier at each step.
# Output:
#' @return Sparsified network, \code{H}, which still maintains evaluator function, \code{func}, plus/minus \code{tol}.
#' @author Alexander Mercier
#' @author Andrew Kramer
#' @examples
#' #Set scoring function
#' mean.weight.degree <- function(graph){
#' graph.ob <- igraph::graph_from_edgelist(graph[,1:2])
#' igraph::E(graph.ob)$weight <- graph[,3]
#' return(mean(igraph::strength(graph.ob)))
#' }
#'
#' #Generate random graph
#' g <- igraph::erdos.renyi.game(100, 0.1)
#' igraph::E(g)$weight <- rexp(length(igraph::E(g)), rate=10) #random edge weights from exp(10)
#' E_List <- cbind(igraph::as_edgelist(g), igraph::E(g)$weight)
#' colnames(E_List) <- c("n1", "n2", "weight")
#' sparse_dist <- simplifyNet::irefit(E_List, func=mean.weight.degree, tol = 0.1)
#' @export
irefit <- function(network, func, tol, rank = 'none', connected = FALSE, directed = FALSE, per = 0.5){
  E_List = simplifyNet::net.as(network, net.to = "E_List", directed = directed)
  if(dim(E_List)[2] == 2){
  } else if(rank == 'weight'){
    E_List = E_List[order(E_List[,3], decreasing = TRUE), ]
  } else if(rank == 'random'){
    E_List = E_List[sample(nrow(E_List)), ]
  } else if(rank == "none"){
    E_List = E_List
  } else {
    E_List = rerank(E_List, rank)
  }
  org_score = func(E_List)
  if(methods::is(org_score, "list")){
    E_List = rerank(E_List, org_score[[2]]) #Refitted weights
    org_score = org_score[[1]]
    stepsize = ceiling(per*nrow(E_List))
    pivot = nrow(E_List) - stepsize
    S_List = E_List[1:pivot,]
    spr_score = func(S_List)
    S_List = rerank(S_List, spr_score[[2]])
    spr_score = spr_score[[1]]
  } else {
    stepsize = ceiling(per*nrow(E_List))
    pivot = nrow(E_List) - stepsize
    S_List = E_List[1:pivot,]
    spr_score = func(S_List)
  }

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

#' @name sparse.step
#' @title Iterative refitting sparsification step
#' @description Recursive sparsification step for iterative refitting.
#' @keywords internal
#' @details Intended as internal function.
#' @return Return sparsified edge list.
#' @author Andrew Kramer
#' @author Alexander Mercier
#' @export
sparse.step <- function(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected){
  S = simplifyNet::EList_Mtrx(S_List, TRUE, n=max(E_List[,1:2]))
  print(sprintf("Stepsize %d, EdgeNum %d, Score %f, Connected %s", stepsize, nrow(S_List), spr_score, as.character(is.connected(S))))
  stepsize = ceiling(per * stepsize)
  if(connected){
    if(abs(spr_score - org_score) > tol || !(is.connected(S)) || spr_score == Inf){
      pivot = nrow(S_List) + stepsize
      S_List = E_List[1:pivot,]
      spr_score = func(S_List)
      if(methods::is(spr_score, "list") && spr_score[[1]] != Inf){
        S_List = rerank(S_List, spr_score[[2]]) #refit
        spr_score = spr_score[[1]] #score
      } else {spr_score = spr_score[[1]]}
      sparse.step(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected)
    } else if(stepsize == 1) {
        return(S_List)
    } else if(abs(spr_score - org_score) <= tol) {
      pivot = nrow(S_List) - stepsize
      S_List = S_List[1:pivot,]
      spr_score = func(S_List)
      if(methods::is(spr_score, "list") && spr_score[[1]] != Inf){
        S_List = rerank(S_List, spr_score[[2]]) #refit
        spr_score = spr_score[[1]] #score
      } else {spr_score = spr_score[[1]]}
      sparse.step(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected)
    }
  } else {
    if(abs(spr_score - org_score) > tol) {
      pivot = nrow(S_List) + stepsize
      S_List = E_List[1:pivot,]
      spr_score = func(S_List)
      if(methods::is(spr_score, "list") && spr_score[[1]] != Inf){
        S_List = rerank(S_List, spr_score[[2]]) #refit
        spr_score = spr_score[[1]] #score
      } else {spr_score = spr_score[[1]]}
      sparse.step(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected)
    } else if(stepsize == 1) {
      return(S_List)
    } else if(abs(spr_score - org_score) <= tol) {
      pivot = nrow(S_List) - stepsize
      S_List = S_List[1:pivot,]
      spr_score = func(S_List)
      if(methods::is(spr_score, "list") && spr_score[[1]] != Inf){
        S_List = rerank(S_List, spr_score[[2]]) #refit
        spr_score = spr_score[[1]] #score
      } else {spr_score = spr_score[[1]]}
      sparse.step(E_List, S_List, stepsize, spr_score, org_score, func, tol, per, connected)
    }
  }
}

#' @name rerank
#' @title Rerank edges
#' @description Rerank edges for iterative refitting.
#' @keywords internal
#' @details Intended as internal function.
#' @return Return an edge list reordered according to \code{rank}.
#' @author Andrew Kramer
#' @author Alexander Mercier
#' @export
rerank <- function(E_List, rank){
  E_List = E_List[order(cbind(E_List, methods::as(rank, 'vector'))[,4], decreasing = TRUE), ]
  return(E_List)
}

#' @name remove.edges
#' @title Remove edges
#' @description remove edges for iterative refitting.
#' @keywords internal
#' @details Intended as internal function.
#' @return Return a truncated edge list with a decreased number of edges determined by \code{per}.
#' @author Andrew Kramer
#' @author Alexander Mercier
#' @export
remove.edges <- function(per, E_List, S_List){
  m = nrow(E_List)
  s = nrow(S_List)
  frac = s/m
  new_frac = ceiling(frac * (1-per) * nrow(E_List))
  S_List = E_List[1:new_frac,]
  return(S_List)
}

#' @name add.edges
#' @title Add edges
#' @description Add edges for iterative refitting.
#' @keywords internal
#' @details Intended as internal function.
#' @return Return an expanded edge list with an increased number of edges determined by \code{per}.
#' @author Andrew Kramer
#' @author Alexander Mercier
#' @export
add.edges <- function(per, E_List, S_List) {
  m = nrow(E_List)
  s = nrow(S_List)
  frac = s/m
  new_frac = ceiling((frac * (1-per) + 1) * nrow(S_List))
  S_List = E_List[1:new_frac,]
  return(S_List)
}

#' @name find.neighbors
#' @title Find Neighbors
#' @description  Finds the neighbors of a given node, \code{v}.
# Input:
#' @param Adj
#' Adjacency matrix
#' @param v
#' Node to find the neighbors of.
# Output:
#' @return Integers designating node indices of the adjacency matrix for the neighbors of \code{v}.
#' @keywords internal
#' @details Intended as internal function.
#' @author Andrew Kramer
#' @author Alexander Mercier
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
#' Adjacency matrix.
#' @param v
#' Node to perform DFS from.
#' @param discovered
#' A list of discovered nodes from \code{v}. If initializing the search, should be all FALSE.
# Output:
#' @return Logical of length n where TRUE denotes connected to node \code{v}.
#' @keywords internal
#' @details Intended as internal function.
#' @author Andrew Kramer
#' @author Alexander Mercier
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
#' Adjacency matrix.
# Output:
#' @return Return TRUE if network is connected and FALSE if not connected.If the network is directed, then this function checks if the network is strongly connected.
#' @keywords internal
#' @details Intended as internal function.
#' @author Andrew Kramer
#' @author Alexander Mercier
#' @export
is.connected <- function(Adj){
  discovered = matrix(FALSE, nrow=1, ncol=ncol(Adj))
  v = sample(1:ncol(Adj), 1)
  discovered = DFS(Adj, v, discovered)
  return(all(discovered))
}




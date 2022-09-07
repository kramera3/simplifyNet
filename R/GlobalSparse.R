#' @name gns
#' @title Global Network Sparsification
#' @description Remove all edges under certain edge weight threshold.
#  Input:
#' @param  network
#' Weighted adjacency matrix, weighted \code{igraph} network, or edge list formatted | n1 | n2 | weight | with colnames \code{c("n1", "n2", "weight")}.
#' @param remove.prop
#' The proportion of highest weighted edges to retain. A value between 0 and 1.
#' @param cutoff
#' Threshold value for edge weight thresholding.
#' @param directed
#' If \code{TRUE}, specifies that the inputted network is directed. Default is \code{FALSE}.
# Output:
#' @return Edge list of sparsified network
#' @author Andrew Kramer
#' @author Alexander Mercier
#' @examples
#' #Generate random ER graph with uniformly random edge weights
#' g = igraph::erdos.renyi.game(100, 0.1)
#' igraph::E(g)$weight <- runif(length(igraph::E(g)))
#' #Sparsify g via GNS
#' S = gns(g, remove.prop = 0.5)
#' sg = simplifyNet::net.as(S, net.to="igraph", directed=FALSE)
#' igraph::ecount(sg)/igraph::ecount(g)#fraction of edges in the sparsifier
#' @export

#Returns the sparse network by some global threshold
gns <- function(network, remove.prop, cutoff, directed = FALSE){
  E_List = simplifyNet::net.as(network, net.to = "E_List", directed = directed)
  if(!missing(remove.prop) && !missing(cutoff)){

    stop("Please only use remove.prop or cutoff.")

  } else if(!missing(remove.prop)){
    E_List = simplifyNet::net.as(network, net.to = "E_List", directed = directed)
    colnames(E_List)<-c("n1","n2","weight")
    ret <- E_List[order(E_List[,3], decreasing=TRUE),]
    ret <- ret[1:ceiling(dim(E_List)[1]*remove.prop),]
  } else if(!missing(cutoff)) {
    E_List = simplifyNet::net.as(network, net.to = "E_List", directed = directed)
    colnames(E_List)<-c("n1","n2","weight")
    ret <- E_List[E_List[,3] %in% E_List[,3][E_List[,3] > cutoff],]
  }

  return(ret)
}

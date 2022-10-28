# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# From Toivonen et al. 2010, "A Framework for Path-Oriented Network Simplification"
# https://link.springer.com/chapter/10.1007/978-3-642-13062-5_21

#' @name bestpath
#' @title Sparsification via Best Path
#' @description Calculates network sparsifier from best path
#  Input:
#' @param  network
#' Weighted adjacency matrix, weighted \code{igraph} network, or edge list formatted | n1 | n2 | weight | with colnames \code{c("n1", "n2", "weight")}.
#' @param directed
#' If \code{TRUE}, specifies that the inputted network is directed. Default is \code{FALSE}.
#' @param associative
#' Designates if the network is associative where edge weight determines "similarity" or "strength" or dissociative where edge weight denotes "dissimilarity" or "distance".\cr
#' If the network is associative, then the shortest path would be found by looking at \code{w_e^-1} where weaker association between nodes suggests a larger distance between nodes for shortest paths. \cr
#' If the network is dissociative, then the shortest path would be between \code{w_e}.
# Output:
#' @return Edge list of sparsified network via best path.
#' @author Alexander Mercier
#' @author Andrew Kramer
#' @references Toivonen, H., Mahler, S., & Zhou, F. (2010, May). A framework for path-oriented network simplification. In International Symposium on Intelligent Data Analysis (pp. 220-231). Springer, Berlin, Heidelberg.
#' @examples
#' #Generate random ER graph with uniformly random edge weights
#' g = igraph::erdos.renyi.game(50, 0.1)
#' igraph::E(g)$weight <- runif(length(igraph::E(g)))
#' #Sparsify g via bestpath
#' S = simplifyNet::bestpath(g, directed = FALSE, associative = TRUE) #Show edge list conversion
#' sg = simplifyNet::net.as(S, net.to="igraph", directed=FALSE)
#' igraph::ecount(sg)/igraph::ecount(g)#fraction of edges in the sparsifier
#' @export

# Returns the sparse network containing every link that is in a best path, Toivonen et al
bestpath <- function(network, directed=FALSE, associative=TRUE){
  ## Determine if network is symmetric or not ##
  grph = simplifyNet::net.as(network, net.to = "igraph", directed = directed)

  if(associative){
    igraph::E(grph)$weight <- (igraph::E(grph)$weight)^(-1) #Change edge weights for best path #Note this is only for wanting the distance for association (weaker association -> larger distance)
  }

  n <- igraph::gorder(grph) # Want number of nodes
  ## Make solution matrix ##
  mask<-matrix(FALSE, nrow=n, ncol=n)
  ## Look for best paths (in each direction) ##
  for (i in 1:n){
    paths <- igraph::get.shortest.paths(grph, from=as.character(i))[[1]] #The initial list has 4 elements, all the best paths are stored in element[[1]] #Does this include weights?
    ## Transform path to "mask" indices of 1 ##
    for(j in 1:length(paths)){
      if(length(paths[[j]]) > 1){ #The path of each to itself is undefined, numeric(0)
        for(k in 1:(length(paths[[j]]) - 1)){
          mask[paths[[j]][k], paths[[j]][k+1]] <- TRUE #From stored edge, e_i - e_{i+1} in the ordered path
        }
      }
    }
  }

  A = as.matrix(simplifyNet::net.as(network, "Adj", directed = directed)) # Return edge list
  ret = mask * A
  ret = simplifyNet::Mtrx_EList(ret, directed = directed)
  return(ret) #Faster way of doing this for large matrices?
}


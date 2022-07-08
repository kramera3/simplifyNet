# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' @author Alexander Mercier
#' @author Andrew Kramer
#' @rdname ToivonenBestPath

# From Toivonen et al. 2010, "A Framework for Path-Oriented Network Simplification"
# https://link.springer.com/chapter/10.1007/978-3-642-13062-5_21

#' @name toivonen
#' @title Sparsification via Best Path
#' @description Calculates network sparsifier from best path
#  Input:
#' @param  E_List
#' Edge list of the given network in the format of | node 1 | node 2 | weight |.
#' @param directed
#' Specifies if the network is directed or undirected. Default is set to undirected.
#' @param associative
#' Designates if the network is associative where edge weight determines "similarity" or "strength" or dissociative where edge weight denotes "dissimilarity" or "distance".\cr
#' If the network is associative, then the shortest path would be found by looking at \code{w_e^-1} where weaker association between nodes suggests a larger distance between nodes for shortest paths.
#' If the network is dissociative, then the shortest path would be between \code{w_e}.
# Output:
#' @return Edge list of sparsified network
#' @examples
#' g = igraph::erdos.renyi.game(100, 0.1)
#' igraph::E(g)$weight <- runif(length(igraph::E(g)))
#' E_List = cbind(igraph::as_edgelist(g), igraph::E(g)$weight)
#' H = toivonen(E_List, directed = FALSE, associative = TRUE)
#' @export

# Returns the sparse network containing every link that is in a best path, Toivonen et al
toivonen <- function(E_List, directed=FALSE, associative=TRUE){
  ## Determine if network is symmetric or not ##

  colnames(E_List) <- c("n1", "n2", "weight")
  df <- data.frame(E_List)
  grph <- igraph::graph.data.frame(df, directed = directed)

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

  E_Mask = simplifyNet::Mtrx_EList(mask, directed = directed) # Return edge list
  ret = merge(E_List, E_Mask, by=c("n1", "n2"))
  ret[,4] <- NULL
  colnames(ret) <- c("n1", "n2", "weight")
  return(ret) #Faster way of doing this for large matrices?
}


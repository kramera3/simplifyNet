# From Foti et al. 2011 "Nonparametric Sparsification of Complex Multiscale Networks"
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016431#s5

#' @name lans
#' @title Local Adaptive Network Sparsification
#' @description Remove all edges under certain probability of the fractional edge weight, \code{alpha}.
#  Input:
#' @param  network
#' Weighted adjacency matrix, weighted \code{igraph} network, or edge list formatted | n1 | n2 | weight | with colnames \code{c("n1", "n2", "weight")}.
#' @param alpha
#' The \code{alpha} value is a predetermined threshold to designate statistically important edges by their fractional edge weight at each node. If the probability of choosing that edge via the CDF is less than or equal to \code{alpha}, then the edge is not included.
#' @param output
#' If the output should be directed or undirected. Default is that the output is the same as the input based on adjacency matrix symmetry. If the default is overridden, set as either "undirected" or "directed".
#' @param directed
#' If \code{TRUE}, specifies that the inputted network is directed. Default is \code{FALSE}.
# Output:
#' @return Weighted adjacency matrix of sparsified network.
#' @author Andrew Kramer
#' @author Alexander Mercier
#' @details For more information on finding alpha values, see: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016431#s5
#' @references Foti, N. J., Hughes, J. M., & Rockmore, D. N. (2011). Nonparametric sparsification of complex multiscale networks. PloS one, 6(2), e16431.
#' @examples
#' #Generate random ER graph with uniformly random edge weights
#' g = igraph::erdos.renyi.game(100, 0.1)
#' igraph::E(g)$weight <- runif(length(igraph::E(g)))
#' #Sparsify g via LANS
#' S = lans(g, alpha = 0.3, output = "undirected", directed = FALSE)
#' #Convert sparsifier to edge list
#' S_List = simplifyNet::Mtrx_EList(S, directed = FALSE)
#' sg = simplifyNet::net.as(S_List, net.to="igraph", directed=FALSE)
#' igraph::ecount(sg)/igraph::ecount(g)#fraction of edges in the sparsifier
#' @export


#Returns the sparse network based on Locally adaptive network sparsification, _F, modified for directed network
lans <- function(network, alpha, output, directed = FALSE){
  Adj = simplifyNet::net.as(network, net.to="Adj", directed=directed)
  #Determine if the outputted network should be directed or not
  if(!missing(output)){
    if(output == "undirected"){
      directed <- FALSE
    } else if(output == "directed"){
      directed <- TRUE
    }
  } else {
    if(Matrix::isSymmetric(Adj)){
      directed <- FALSE
    } else {
      directed <- TRUE
    }
  }
  n <- dim(Adj)[1]
  mask<-matrix(0L, nrow=n, ncol=n)
  Adj = as.matrix(Adj)
  ##Get the fractional weights of the edges for each node
  ##go across rows and divide by sum, then across columns, calculating one triangle at a time
  ##In the models source is j and destination is i, that explains switch of rows and columns with Foti et al Supplement S4
  Pij<-proportions(Adj, margin=1) #Relative prob of numeric in the matrix across rows.
  lower<-proportions(Adj, margin=2) #Relative prob of numeric in the matrix across columns. #Note# This will be the same as above iff the matrix is symmetric
  Pij[lower.tri(Pij)]<-lower[lower.tri(lower)]# Make sure correct prob is in lower triangle
  if (length(alpha>0)){
    if(!directed){ #If ouputted network is undirected
      for (i in 1:n){
        temp<-stats::ecdf(Pij[,i]) #apply or for each parallel
        pa<-stats::quantile(temp,1-alpha)
        for (j in 1:n){
          if (Pij[j,i]>pa) {
            mask[j,i]<-1 #Creates a mask
            mask[i,j]<-1 #Undirected mask
          }
        }
      }
    } else { #If outputted network is directed
      for (i in 1:n){
        temp<-stats::ecdf(Pij[,i]) #apply or for each parallel
        pa<-stats::quantile(temp,1-alpha)
        for (j in 1:n){
          if (Pij[j,i]>pa) {
            mask[j,i]<-1 #Creates a mask
          }
        }
      }
    }
  }

  return(mask * Adj)

}

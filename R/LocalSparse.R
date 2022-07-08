#' @author Andrew Kramer
#' @author Alexander Mercier
#' @rdname LocalSparse

# From Foti et al. 2011 "Nonparametric Sparsification of Complex Multiscale Networks"
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016431#s5

#' @name lans
#' @title Global threshold sparsification
#' @description Remove all edges under certain edge weight threshold.
#' Used to fulfill Johnson-Lindenstrauss lemma \cr
#' From Spielman and Srivastava 2011 \cr
#  Input:
#' @param  Adj
#' Weighted adjacency matrix of the network of class "matrix".
#' @param remove.prop
#' Alpha value threshold to designate statistically "unimportant" edges by edge weight.
#' @param output
#' If the output should be directed or undirected. Default is that the output is the same as the input based on adjacency matrix symmetry. If the default is overridden, set as either "undirected" or "directed".
# Output:
#' @return Weighted adjacency matrix of sparsified network
#' @examples
#' g = igraph::erdos.renyi.game(100, 0.1)
#' igraph::E(g)$weight <- runif(length(igraph::E(g)))
#' A = igraph::as_adj(g, attr = "weight")
#' A = as(A, "matrix")
#' H = lans(A, remove.prop = 0.3)
#' @export


#Returns the sparse network based on Locally adaptive network sparsification, _F, modified for directed network
lans <- function(Adj, remove.prop, output){

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
  ##Get the fractional weights of the edges for each node
  ##go across rows and divide by sum, then across columns, calculating one triangle at a time
  ##In the models source is j and destination is i, that explains switch of rows and columns with Foti et al Supplement S4
  Pij<-proportions(Adj, margin=1) #Relative prob of numeric in the matrix across rows.
  lower<-proportions(Adj, margin=2) #Relative prob of numeric in the matrix across columns. #Note# This will be the same as above iff the matrix is symmetric
  Pij[lower.tri(Pij)]<-lower[lower.tri(lower)]# Make sure correct prob is in lower triangle
  if (length(remove.prop>0)){
    if(!directed){ #If ouputted network is undirected
      for (i in 1:n){
        temp<-stats::ecdf(Pij[,i]) #apply or for each parallel
        pa<-stats::quantile(temp,1-remove.prop)
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
        pa<-stats::quantile(temp,1-remove.prop)
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

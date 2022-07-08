#' @author Alexander Mercier
#' @rdname GlobalSparse


#' @name gns
#' @title Global Threshold Sparsification
#' @description Remove all edges under certain edge weight threshold.
#  Input:
#' @param  E_List
#' Edge list of the given network in the format of | node 1 | node 2 | weight |.
#' @param remove.prop
#' The proportion of highest weighted edges to retain. A value between 0 and 1.
#' @param cutoff
#' Threshold value for edge weight thresholding.
# Output:
#' @return Edge list of sparsified network
#' @examples
#' g = igraph::erdos.renyi.game(100, 0.1)
#' igraph::E(g)$weight <- runif(length(igraph::E(g)))
#' E_List = cbind(igraph::as_edgelist(g), igraph::E(g)$weight)
#' H = gns(E_List, remove.prop = 0.5)
#' @export

#Returns the sparse network by some global threshold
gns <- function(E_List, remove.prop, cutoff){

  if(!missing(remove.prop) && !missing(cutoff)){

    stop("Please only use remove.prop or cutoff.")

  } else if(!missing(remove.prop)){

    colnames(E_List)<-c("n1","n2","weight")
    ret <- E_List[order(E_List[,3], decreasing=TRUE),]
    ret <- ret[1:ceiling(dim(E_List)[1]*remove.prop),]

  } else if(!missing(cutoff)) {

    colnames(E_List)<-c("n1","n2","weight")
    ret <- E_List[E_List[,3] %in% E_List[,3][E_List[,3] > cutoff],]

  }

  return(ret)
}

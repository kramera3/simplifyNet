#Input: Adj matrix or Edge List (| from | to | weight |)
#Output: Edge List (| from | to | weight |)


# Network sparsification
# Possible methods
# 1. Toivonen - Remove all edges that are not part of best path, algorithm from Toivonen et al. implemented with 'igraph' package
# 2. GlobalSparse - Remove weakest edges globally, default is to refit after each until reaching statistical cutoff, can also set proportion or go until disconnected
# 3. LocalAdapt - Remove weakest edges locally rather than globally, can set proportion or use refitting and cutoff

simplifyNet <- function(data, func, method="Toivonen", weights=NULL, cutoff, remove.prop){
  # Removes edges from network models using algorithm of choice
  #
  # Args:
  #   model: function that can fit the network of interest, should generally return vector of fit parameters or matrix of edge weights, as well negative log-likelihood or some measure to use for cutoff
  #   data: list of named data structures needed by 'model', need to match objects called in model
  #         For gravity model this includes distances, mass, for White-nose syndrome also need winter length, species richness
  #   parameters: fit parameters for the full (initial) network model if model is not meant to be fit on the fly
  #   method: algorithm used to remove edges, currently includes Toivonen (best-path), Refit and LocalAdapt. Default is Toivonen.
  #   weights: directly pass the edge weights, can save steps in process but only sufficient for Toivonen
  #   cutoff: statistical cutoff used to halt edge removal, default is 2 (NLL units)
  #   search: whether to remove chunks of edges or proceed one-at-a-time. "direct" refits after each edge, can often be computationally prohibitive
  #   disconnect: whether to allow sparsification to lead to disconnection, only relevant for GlobalSparse
  #   remove.prop: argument can be passed to LocalAdapt or GlobalSparse to remove a set proportion of edges instead of using statistical cutoff

  if (method == "Toivonen"){
    output <- toivonen(data=data)
    output <- adj_to_edgelst(output)
    edge_list <- output
  }
  if (method == "GlobalSparse"){
    output <- gns(data=data, cutoff=cutoff)
    edge_list <- output
  }
  if (method == "LocalAdapt"){
    output <- lans(data=data, remove.prop=remove.prop)
    output <- adj_to_edgelst(output)
    edge_list <- output
  }
  return(edge_list)
}

# Returns the sparse network containing every link that is in a best path, Toivonen et al
toivonen <- function(data){
  ## Make sure edgelist is in adj matrix format ##
  if(is.matrix(data)!= TRUE) {
    colnames(data) <- c("from", "to", "weight")
    data <- get.adjacency(graph.data.frame(data), attr = "weight", sparse = FALSE)
  }
  ## Create graph from adj matrix for get shortest path ##
  grph <- graph_from_adjacency_matrix(data, mode = "directed", weighted = TRUE)

  n <- dim(data)[1]
  ## Make solution matrix ##
  mask<-matrix(FALSE, nrow=n, ncol=n)
  ## Look for best paths (in each direction) ##
  for (i in 1:n){
    paths <- get.shortest.paths(grph, from=as.character(i))[[1]] #The initial list has 4 elements, all the best paths are stored in element[[1]]
    ## Transform path to "mask" indices of 1 ##
    for(j in 1:length(paths)){
      if(length(paths[[j]]) > 1){ #The path of each to itself is undefined, numeric(0)
        for(k in 1:(length(paths[[j]]) - 1)){
          mask[paths[[j]][k], paths[[j]][k+1]] <- TRUE
        }
      }
    }
  }
  return(mask * data)
}

#Returns the sparse network by some global threshold
gns <- function(data, cutoff){
  if (is.matrix(data) == TRUE) {
    data <- adj_to_edgelst(data)
  }
  if (cutoff > 0) {
    colnames(data)<-c("from","to","weight")
    ret <- data[data$weight %in% data$weight[data$weight > cutoff],]
    return(ret)
  }
}

#Returns the sparse network based on Locally adaptive network sparsification, Foti et al 2011, modified for directed network
lans <- function(data, remove.prop){
  if(is.matrix(data)!= TRUE) {
    colnames(data) <- c("from", "to", "weight")
    data <- get.adjacency(graph.data.frame(data), attr = "weight", sparse = FALSE)
  }
  n <- dim(data)[1]
  Aij<-matrix(0, nrow=n, ncol=n)
  ##Get the fractional weights of the edges for each node
  ##go across rows and divide by sum, then across columns, calculating one triangle at a time
  ##In the models source is j and destination is i, that explains switch of rows and columns with Foti et al Supplement S4
  Pij<-prop.table(data,margin=1) #Relative prob of numeric in the matrix across rows.
  lower<-prop.table(data, margin=2) #Relative prob of numeric in the matrix across columns. #Note# This will be the same as above iff the matrix is symetric
  Pij[lower.tri(Pij)]<-lower[lower.tri(lower)]#lower.tri gives a matrix of same size with TRUE in the lower triangle (defult is not to include diagonal)
  if (length(remove.prop>0)){
    for (i in 1:n){
      temp<-ecdf(Pij[,i])
      pa<-quantile(temp,1-remove.prop) #?#
      for (j in 1:n){
        if (Pij[j,i]>pa) Aij[j,i]<-1
      }
    }
  } else { # call the iterative refitting algorithm that checks against cutoff
    Aij<-Pij # Placeholder #?# For what? ^
  }
  return(Aij*data)
}

adj_to_edgelst<-function(mtrx,weights=TRUE){
  adj <- graph.adjacency(mtrx,weighted=weights)
  ej_lst <- get.data.frame(adj)
  return(ej_lst)
}


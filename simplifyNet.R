# Author: Andrew M. Kramer

# Script for network sparsification
# Will have a single wrapper function that calls options for each version
# Will rely on supplied functions from outside for refitting, returning negative log-likelihood (or similar), these functions should be sourced and available in the environment
# Inputs are data for fitting full network, parameters estimated for full network
# User also supplies some options

# Network sparsification
# Possible methods
  # 1. Toivonen - Remove all edges that are not part of best path, algorithm from Toivonen et al. implemented with 'igraph' package
  # 2. GlobalSparse - Remove weakest edges globally, default is to refit after each until reaching statistical cutoff, can also set proportion or go until disconnected
  # 3. LocalAdapt - Remove weakest edges locally rather than globally, can set proportion or use refitting and cutoff

simplifyNet <- function(data, parameters=NULL, method="Toivonen", weights=NULL, cutoff=0.1, search="adaptive", disconnect=F, remove.prop=NULL){
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
    output <- toivonen(data=data, parameters=parameters, weights=weights)
    edges <- output
  }
  if (method == "GlobalSparse"){
    output <- gns(data=data, parameters=parameters, weights=weights, cutoff=cutoff, search=search, disconnect=disconnect, remove.prop=remove.prop)
    edges <- output[[1]]
    parameters <- output [[2]]
  }
  if (method == "LocalAdapt"){
    output <- lans(data=data, parameters=parameters, weights=weights, cutoff=cutoff, remove.prop=remove.prop)
    edges <- output [[1]]
    parameters <- output [[2]]
  }
}

# Returns the sparse network containing every link that is in a best path, Toivonen et al 
toivonen <- function(data, parameters, weights){
  beta <- parameters
  require(igraph)
  if (length(weights) == 0){
    weights <- getweights(dist=data$dist,mass=data$mass,tau=data$tau,SR=data$SR,beta) ##
  }
  weights <- as.matrix(weights)
  n <- dim(weights)[1]
  nlweights <- -log(weights)
  colnames(nlweights) <- as.character(1:n)
  graph <- graph.adjacency(nlweights, mode="directed", weighted=TRUE)
  
  ###Make solution matrix
  mask<-matrix(FALSE, nrow=n, ncol=n)
  
  ###Look for best paths (in each direction)
  for (i in 1:n){
    paths <- get.shortest.paths(graph, from=as.character(i))[[1]] ##The initial list has 4 elements, all the best paths are stored in element[[1]]
    ##transform path to "mask" indices of 1
    for(j in 1:length(paths)){
      if(length(paths[[j]]) > 1){ ##The path of each to itself is undefined, numeric(0)
        for(k in 1:(length(paths[[j]]) - 1)){
          mask[paths[[j]][k], paths[[j]][k+1]] <- TRUE
        }
      }
    }
  }
  return(mask)
}

# Returns the sparse network estimated through dropping weak links and refitting. 
gns <- function(data, parameters, cutoff, search, disconnect, remove.prop){
  beta <- parameters
  dist<-data$dist
  if (length(weights) == 0){
    weights <- getweights(dist,data$mass,data$tau,data$SR,beta)
  }
  weights <- as.matrix(weights)
  if (length(remove.prop)>0){
    d <- dist
    mask <- d!=0
    w <- weights*mask
    o <- order(w)
    d[o[1:floor(remove.prop*length(o))]] <- 0
    mask <- d!=0
    if (edge_connectivity(mask) == 0) return (print("Resulting network is not strongly connected. Set 'disconnect=F' to ignore"))
  } else {
    while(1){  
      i=sum(dist==0)
      print(i)
      dist=rmbinary(dist, data$mass, data$tau, data$SR, beta, data, years) #Years?
      tbeta=try(refit2(dist,mass,tau,SR,data,beta,years),silent=TRUE)
      if(is.character(tbeta)){break}
      beta=tbeta
      mask=dist!=0
      save(mask, beta, file=paste('relaxbinary', i, '.Rdata', sep=''))
    }
  }
}

#Returns the sparse network based on Locally adaptive network sparsification, Foti et al 2011, modified for directed network
lans <- function(data, parameters, weights, cutoff, remove.prop){
  beta <- parameters
  if (length(weights) == 0){
    weights <- getweights(data$dist,data$mass,data$tau,data$SR,beta)
  }
  weights <- as.matrix(weights)
  n <- dim(weights)[1]
  Aij<-matrix(0, nrow=n, ncol=n)
  ##Get the fractional weights of the edges for each node
  ##go across rows and divide by sum, then across columns, calculating one triangle at a time
  ##In the models source is j and destination is i, that explains switch of rows and columns with Foti et al Supplement S4
  Pij<-prop.table(weights,margin=1)
  lower<-prop.table(weights, margin=2)
  Pij[lower.tri(Pij)]<-lower[lower.tri(lower)]
  if (length(remove.prop>0)){
    for (i in 1:n){
      temp<-ecdf(Pij[,i])
      pa<-quantile(temp,1-remove.prop)
      for (j in 1:n){
        if (Pij[j,i]>pa) Aij[j,i]<-1
      }
    }
  } else { # call the iterative refitting algorithm that checks against cutoff
    Aij<-Pij # Placeholder
  }
  return(Aij)
}



#partial relaxation with binary search
relaxbinary<-function(dist, mass, tau, SR, beta, data, years){
  while(1){  
    i=sum(dist==0)
    print(i)
    dist=rmbinary(dist, mass, tau, SR, beta, data, years)
    tbeta=try(refit2(dist,mass,tau,SR,data,beta,years),silent=TRUE)
    if(is.character(tbeta)){break}
    beta=tbeta
    mask=dist!=0
    save(mask, beta, file=paste('relaxbinary', i, '.Rdata', sep=''))
  }
}

#use binary search to find nll threshold faster
rmbinary<-function(dist, mass, tau, SR,beta, data, years, remove.prop){
  dist=as.matrix(dist)
  mask=dist!=0
  w=getweights(dist, mass, tau, SR,beta)
  w=w*mask
  o=order(w)
  beginll=objf(beta, dist, mass, tau, SR, data, years)
  lower=sum(dist==0)+1
  upper=prod(dim(dist))
  if(lower>upper){lower=upper}
  while(lower != upper){
    d=dist
    pivot=floor(mean(c(upper,lower)))
    d[o[1:pivot]]<-0
    testnll=signif(beginll,7)==signif(objf(beta, d, mass, tau, SR, data, years), 7)
    if(testnll){lower=pivot+1}
    else{upper=pivot}
  }
  d=dist
  d[o[1:pivot]]<-0
  return(d)
}

#Spatial decay function
decay<-function(x,mass,tau,SR,beta0,beta1,beta2,beta3,beta4){
  return(tmp=1/(1+exp(beta0 +beta1*tau + (beta2*x)/(mass**beta3) + SR*beta4)))
}

getweights<-function(dist,mass,tau,SR,beta){
  mass=mass%*%t(mass)
  tau=matrix(tau,ncol=length(tau),nrow=length(tau),byrow=TRUE)
  SR=matrix(SR,ncol=length(SR),nrow=length(SR),byrow=TRUE)
  ret=decay(dist,mass,tau,SR,beta[1],beta[2],beta[3],beta[4],beta[5])
  diag(ret)<-0
  return(ret)
}

refit2<-function(dist,mass,tau,SR,data,start,years){
  last<-c(0,0,0,0)
  beta<-start
  while(!all(signif(beta,7)==signif(last,7))){
    last<-beta
    beta<-optim(par=last,fn=objf,dist=dist,mass=mass,tau=tau,SR=SR,data=data,years=years)$par
  }
  return(beta)
}

#returns probabilities that each node will be infected in the next time step
getprobs<-function(start,dist,mass,tau,SR,b){
  
  #Check Trivial Cases
  if(all(start)){
    return(rep(1,length(start)))
  }
  if(!any(start)){
    return(rep(0,length(start)))
  }
  
  #calculate probability of infection from each infected node node to each uninfected node
  mass=mass%*%t(mass)
  tmp=decay(dist[which(start),which(!start)],mass[which(start),which(!start)], (rep(1,length(which(start))) %*% t(tau[which(!start)])), (rep(1,length(which(start))) %*% t(SR[which(!start)])), b[1],b[2],b[3],b[4],b[5])
  tmp=tmp*(dist[which(start),which(!start)] != 0)
  
  #combine probabilities of infected nodes for each uninfected node
  if(!is.null(dim(tmp))){
    tmp=log(1-tmp)
    tmp=rep(1,dim(tmp)[1]) %*% tmp
    tmp=1-exp(tmp[1,])
  }else if(sum(!start)==1){
    tmp=1-exp(sum(log(1-tmp)))
  }
  
  #construct return vector
  ret=vector(length=length(start))
  ret[which(start)]=1
  ret[which(!start)]=tmp
  
  return(ret)
}

#calculate log likelihood for a single timestep
lglike<-function(p,d){
  for(i in 1:length(d)){
    if(!d[i]){
      p[i]=1-p[i]
    }
  }
  return(sum(log(p)))
}


##This does the LANS method but check for statistical equivalence as in the relax binary method above
relaxLANS<-function(dist, mass, tau, SR, beta, alpha=0.25){
  while(1){  
    i=sum(dist==0)
    print(i)
    dist=rmbinary(dist, mass, tau, SR, beta, data, years)
    tbeta=try(refit2(dist,mass,tau,SR,data,beta,years),silent=TRUE)
    if(is.character(tbeta)){break}
    beta=tbeta
    mask=dist!=0
    save(mask, beta, file=paste('LANS', i, '.Rdata', sep=''))
  }
}



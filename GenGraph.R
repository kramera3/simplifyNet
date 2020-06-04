gengraph <- function(node, m.log, sd.log, rate, alpha, beta, u.bound = 1, dstrb.type, directed = FALSE, rho = 1)
{
  #Find number of edges
  num.edge <- (node*(node-1))/2
  if ((dstrb.type == "lognormal")|(dstrb.type == "Lognormal")) {
    #Generate random lognormal distribution
    dstrb <- rlnorm(num.edge, meanlog=m.log, sdlog=sd.log)
  } else if ((dstrb.type == "exponetial")|(dstrb.type == "Exponetial")){
    #Generate random exponetial distribution
    dstrb <- rexp(num.edge, rate)
  } else if ((dstrb.type == "beta")|(dstrb.type == "Beta")) {
    #Generate random exponetial distribution
    dstrb <- rbeta(num.edge, alpha, beta)
  }
  dstrb.norm <- (dstrb/(max(dstrb)))*u.bound
  #Create matrix, 0 diagonal along the middle (fully connected)
  adj <- matrix(data = 0L, nrow = node, ncol = node)
  if (directed == FALSE) { #undirected -- symmetric adj matrix... i.e. a_ij = a_ji
    tick.1 <- 1
    tick.2 <- 1
    for(i in 2:node){
      for(j in 1:tick.1){
        adj[i,j]<- dstrb.norm[tick.2] 
        tick.2 <- tick.2 + 1
      }
      tick.1 <- tick.1 + 1
    }
  adj <- adj + t(adj)  
  } else { #Asymmentric adj matrix (directed) i.e. if a_ij = #, then a_ji = 0
    il<-ceiling(runif(num.edge, min = 1, max = 100))
    jl<-ceiling(runif(num.edge, min = 1, max = 100))
    for(i in 1:(num.edge * rho)){
      adj[il[i],jl[i]]<-dstrb.norm[i]
    }
    
  }
  return(adj) 
}
# irefit WNS test script
library(vegan)
library(Matrix)

uc<-read.csv("C:/Users/henry/OneDrive/Desktop/SimplifyNet2022/WNS/us_data3.csv", header=TRUE) #This contains all necessary climate and simulation data added by Mercier to us_data3_D.csv, checked by A. Kramer
uc$X<-NULL
#Subset county data to only Counties with Caves (cc)
cc=uc[which(uc$caves>0),]
#Get coordinates of county centroids
c=rbind(cc$x,cc$y)/1000  #divide by 1000 to convert to km
#Calculate distance matrix
dist=makedist(c)
#Set betas
beta <- c(6.48503691,-0.01255,0.01317286,0.19064664,0.10679775)

WNSWrapper <- function(E_List){
  uc = read.csv("C:/Users/henry/OneDrive/Desktop/SimplifyNet2022/WNS/us_data3.csv", header=TRUE)
  years=max(cc$MaherWNS[which(cc$MaherWNS!=Inf)])
  uc$X<-NULL
  cc=uc[which(uc$caves>0),]
  beta <- c(6.48503691,-0.01255,0.01317286,0.19064664,0.10679775)
  dist = simplifyNet::EList_Mtrx(E_List, TRUE) #Not efficient, but okay - coerce to symmetric matrix
  dist = as.matrix(dist)
  mask=dist!=0
  w=getweights(dist, cc$caves, cc$tau, cc$sr, beta)
  w=w*mask
  w = simplifyNet::Mtrx_EList(w, TRUE)[,3]
  nll=objf(beta, dist, cc$caves, cc$tau, cc$sr, cc$MaherWNS, years)
  return(list(nll, w))
}

E_List = simplifyNet::Mtrx_EList(dist, directed = TRUE)
sparse_dist = simplifyNet::irefit(E_List, func = WNSWrapper, tol = 2) #generate sparsified dist matrix

# Then run all the WNS items... #

######### WNS fitting #########

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

makedist<-function(coords, cutoff=Inf){
  dmat=sqrt(((coords[1,] %*% t(rep(1,length(coords[1,])))) - (rep(1,length(coords[1,])) %*% t(coords[1,])))**2 + ((coords[2,] %*% t(rep(1,length(coords[2,])))) - (rep(1,length(coords[2,])) %*% t(coords[2,])))**2)

  dmat=dmat*(dmat<cutoff)

  return(Matrix(dmat))
}


objf<-function(beta, dist, mass, tau, SR, data, years){
  ret=0

  for(i in 1:(years)){
    p=getprobs((data<i),dist,mass,tau,SR,beta)
    ret=ret+lglike(p,(data<(i+1)))
  }

  return(-ret)
}

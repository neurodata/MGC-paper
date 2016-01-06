#source("verifyNeighbors.R")

verifyNeighbors <-function(p,thres){
  if (missing(thres)){
    thres=0; 
  }
  n=nrow(p);
  pmin=min(p);
  
  neighbor=which((p-pmin)<=thres);
  #neighbor=matrix(0,length(nind),2);
  #nind=nind[length(nind)];
  
  #   ind=((p-pmin)<thres);
  #   nind=which(ind==TRUE);
  #   nind=nind[length(nind)];
  #   
  #neighbor[,2]=ceiling(nind/n);
  #neighbor[,1]=nind-(neighbor[,2]-1)*(n);
  return(neighbor);
}
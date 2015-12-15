# source("disToRanks.R")

disToRanks <- function(dis) {
  # Transform from distance to ranking, order from 0,...,n-1.
  # For ties, the minimum ranking is used.
  n=nrow(dis);
  disRank=matrix(0,n,n);
  for (i in (1:n)){
    v=dis[,i];
    disRank[,i]=rank(v,ties.method="min")-1;
    #disRank[i,i]=0;
  }
  return(disRank);
}
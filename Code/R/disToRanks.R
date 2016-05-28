# source("disToRanks.R")

disToRanks <- function(dis) {
  # An auxiliary function of LocalCorr that calculates ranking within each column for a distance matrix, order from 1,...,n.
  # For ties, the minimum ranking is used.
  n=nrow(dis);
  disRank=matrix(0,n,n);
  for (i in (1:n)){
    v=dis[,i];
    tmp=rank(v,ties.method="min");
    tu=unique(tmp);
    if (length(tu)!=max(tmp)){
      tu=sort(tu);
      for (j in (1:length(tu))){
        kk=which(tmp==tu[j]);
        tmp[kk]=j;
      }
    }
    disRank[,i]=tmp;
  }
  return(disRank);
}
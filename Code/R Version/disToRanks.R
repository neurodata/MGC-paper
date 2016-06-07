# source("disToRanks.R")

disToRanks <- function(dis) {
  # An auxiliary function that sorts the entries within each column by ascending order.
  #
  # The input is assumed to be a distance matrix.
  #
  # The output is column-wise rank, ordered from 1,...,n.
  #
  # For ties, the minimum ranking is used, 
  # e.g. if there are repeating distance entries, the order is like 1,2,3,3,4,..,n-1.
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
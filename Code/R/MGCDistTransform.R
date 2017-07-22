MGCDistTransform <- function(X,Y,option){
  
  if (missing(option)){
    option='mcor'; # use mcorr by default
  }

  # depending on the choice of the global correlation, properly center the distance matrices
  tA=DistCentering(X,option);
  A=tA$A;
  RX=tA$RX;
  tA=DistCentering(Y,option);
  B=tA$A;
  RY=tA$RX;
  
  if (option!='mcorDouble' && option!='dcorDouble'){
    B=t(B);
    RY=t(RY);
  } 
  
  result=list(A=A,RX=RX,B=B,RY=RY);
  return(result);
}

DistCentering<-function(X,option){
  # An auxiliary function that properly centers the distance matrix X,
  # depending on the choice of global corr.
  n=nrow(X);
  RX=DistRanks(X);
  
  # centering for distance correlation / modified distance correlation /
  # Mantel coefficient
  if (option=='dcor'){ # single centering of dcor
    EX=t(matrix(rep(colMeans(X),n), ncol = n));
  }
  if (option=='mcor'){ # single centering of mcor
    EX=t(matrix(rep(colMeans(X),n), ncol = n));
    EX=EX+X/n;
  }
  if (option=='mgc'){ 
    EX=t(matrix(rep(colMeans(X)*n/(n-1),n), ncol = n));
  }
  if (option=='mantel'){
    EX=sum(X)/n/(n-1);
  }
  if (option=='dcorDouble'){ # original double centering of dcor
    EX=t(matrix(rep(colMeans(X),n), ncol = n))+matrix(rep(rowMeans(X),n), ncol = n)-mean(X);
  }
  if (option=='mcorDouble'){ # original double centering of mcor
    EX=t(matrix(rep(colMeans(X),n), ncol = n))+matrix(rep(rowMeans(X),n), ncol = n)-mean(X);
    EX=EX+X/n;
  }
  A=X-EX;
  
  # for mcor or Mantel, exclude the diagonal entries
  # This is a simpler diagonal modification than the original mcor
  if (option!='mcorDouble' && option!='dcorDouble'){
    for (j in (1:n)){
      A[j,j]=0;
    }
  }
  result=list(A=A,RX=RX);
  return(result);
}

DistRanks <- function(dis) {
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
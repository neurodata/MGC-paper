DistCentering<-function(X,option){
  # An auxiliary function that properly centers the distance matrix X,
  # depending on the choice of global corr.
  n=nrow(X);
  
  # centering for distance correlation / modified distance correlation /
  # Mantel coefficient
  if (option=='dcor'){ # single centering of dcor
    EX=t(matrix(rep(colMeans(X),n), ncol = n));
  }
  if (option=='mcor'){ # single centering of mcor
    EX=t(matrix(rep(colMeans(X),n), ncol = n));
    EX=EX+X/n;
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
  if (option=='mcor'||option=='mantel'){
    for (j in (1:n)){
      A[j,j]=0;
    }
  }
  return(A);
}
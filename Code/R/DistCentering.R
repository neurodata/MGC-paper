DistCentering<-function(X,option){
  # An auxiliary function that properly centers the distance matrix X,
  # depending on the choice of global corr.
  n=nrow(X);
  
  # centering for distance correlation / modified distance correlation /
  # Mantel coefficient
  if (option=='dcor'){
    EX=t(matrix(rep(colMeans(X),n), ncol = n));
  }
  if (option=='mcor'){
    EX=t(matrix(rep(colMeans(X),n), ncol = n))+X/n;
  }
  if (option=='mantel'){
    EX=sum(X)/n/(n-1);
  }
  A=X-EX;
  
  # for mcorr or Mantel, exclude the diagonal entries
  if (option=='mcor'||option=='mantel'){
    for (j in (1:n)){
      A[j,j]=0;
    }
  }
  return(A);
}
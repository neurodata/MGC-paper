source("disToRanks.R")

LocalCorr <- function(X,Y,option, disRank){
  # Author: Cencheng Shen
  # Implement the local correlation coefficients from Shen, Jovo, CEP 2016.
  #
  # By specifying option=1, 2, or 3, it calculates the local correlations of dcorr, mcorr, and Mantel
  # Pre-specify a size n * 2n rank matrix can save the sorting for X and Y
  if (missing(option)){
    option=1; # by default use dcorr
  }
  if (missing(disRank)){
    disRank=cbind(disToRanks(X), disToRanks(Y)); # Sort distances within columns, if the ranks are not provided
  }
  n=nrow(X);
  
  # depending on the choice of the global test, calculate the entries of A and B
  # accordingly for late multiplication.
  A=Centering(X,option);
  B=Centering(Y,option);
  
  RX=disRank[,1:n]; # the ranks for X
  RY=disRank[,(n+1):(2*n)]; # the ranks for Y
  result=LocalComputation(t(A),B,RX,RY);
  return(result);
}

Centering<-function(X,option){
  # An auxiliary function that centers the distance matrix X, for dcorr / mcorr / Mantel
  n=nrow(X);
  if (option==3){
    # centering for Mantel
    EX=sum(X)/n/(n-1);
    A=X-EX;
    # Mantel does not use diagonal entries, which is equivalent to set them zero
    for (j in (1:n)){
      A[j,j]=0;
    }
  } else {
    # centering for mcorr/dcorr
    A=X-t(matrix(rep(colMeans(X),n), ncol = n));
    # for mcorr, further adjust the centered matrices to remove high-dimensional bias
    if (option==2){
      A=A-X/n;
      # the diagonals of mcorr are set to zero, instead of the original formulation of mcorr      for (j in (1:n)){
      for (j in (1:n)){
        A[j,j]=0;
      }
    }
  }
  return(A);
}

LocalComputation <- function(A,B,RX,RY){
  # An auxiliary function that computes all local correlations simultaneously
  n=nrow(A);
  nX=max(RX);
  nY=max(RY);
  corrXY=matrix(0,nX,nY);
  varX=rep(0,nX);
  varY=rep(0,nY);
  EX=rep(0,nX);
  EY=rep(0,nY);
  
  # summing up the entriwise product of A and B based on the ranks, which
  # yields the local family of covariance and variances
  for (j in (1:n)){
    for (i in (1:n)){
      a=A[i,j];
      b=B[i,j];
      k=RX[i,j];
      l=RY[i,j];
      corrXY[k,l]=corrXY[k, l]+a*b;
      varX[k]=varX[k]+a^2;
      varY[l]=varY[l]+b^2;
      EX[k]=EX[k]+a;
      EY[l]=EY[l]+b;
    }
  }
  
  for (k in (1:(nX-1))){
    corrXY[k+1,1]=corrXY[k,1]+corrXY[k+1,1];
    varX[k+1]=varX[k]+varX[k+1];
    EX[k+1]=EX[k]+EX[k+1];
  }
  for (l in (1:(nY-1))){
    corrXY[1,l+1]=corrXY[1,l]+corrXY[1,l+1];
    varY[l+1]=varY[l]+varY[l+1];
    EY[l+1]=EY[l]+EY[l+1];
  }
  for (l in (1:(nY-1))){
    for (k in (1:(nX-1))){
      corrXY[k+1,l+1]=corrXY[k+1,l]+corrXY[k,l+1]+corrXY[k+1,l+1]-corrXY[k,l];
    }
  }
  
  # normalize the covariance by the variances yields the local family of correlation.
  options(warn=-1);
  corrXY=(corrXY-EX%*%t(EY)/n^2);
  varX=varX-EX^2/n^2;
  varY=varY-EY^2/n^2;
  corrXY=corrXY/Re(sqrt(varX%*%t(varY))); 
  options(warn=0);
  
  # set any local correlation to 0 if any corresponding local variance is no larger than 0
  for (k in (1:nX)){
    if (varX[k]<=0){
      corrXY[k,]=rep(0,nY);
    }
  }
  for (l in (1:nY)){
    if (varY[l]<=0){
      corrXY[,l]=rep(0,nX);
    }
  }
  result=list(corr=corrXY,varX=varX,varY=varY);
  return(result);
}
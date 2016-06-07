source("disToRanks.R")

LocalCorr <- function(X,Y,option){
  # Author: Cencheng Shen
  # The main function that calculates all local correlation coefficients.
  #
  # The inputs are: 
  # two distance matrices X and Y;
  # an option that specifies which global correlation to use, set to 1,2,3 for dcorr / mcorr / Mantel.
  #
  # The outputs are all local correlations and all local variances.
  if (missing(option)){
    option=2; # use mcorr by default
  }
  n=nrow(X);
  disRank=cbind(disToRanks(X), disToRanks(Y)); # sort distances within columns
  
  # depending on the choice of the global correlation, properly center the distance matrices
  A=Centering(X,option);
  B=Centering(Y,option);
  
  RX=disRank[,1:n]; # the ranks for X
  RY=disRank[,(n+1):(2*n)]; # the ranks for Y
  result=LocalComputation(t(A),B,RX,RY); # compute all local corr / var statistics
  return(result);
}

Centering<-function(X,option){
  # An auxiliary function that properly centers the distance matrix X,
  # depending on the choice of global corr.
  n=nrow(X);
  if (option==3){
    # centering for Mantel
    EX=sum(X)/n/(n-1);
    A=X-EX;
    for (j in (1:n)){
      A[j,j]=0;
    }
  } else {
    # centering for dcorr / mcorr, which uses single centering rather than
    # the original double centering
    A=X-t(matrix(rep(colMeans(X),n), ncol = n));
    
    # for mcorr, further adjust the centered matrices to remove high-dimensional bias
    if (option==2){
      A=A-X/n;
      for (j in (1:n)){
        A[j,j]=0;
      }
    }
  }
  return(A);
}

LocalComputation <- function(A,B,RX,RY){
  # An auxiliary function that computes all local correlations simultaneously in O(n^2)
  n=nrow(A);nX=max(RX);nY=max(RY);
  corrXY=matrix(0,nX,nY);varX=rep(0,nX);varY=rep(0,nY);
  EX=rep(0,nX);EY=rep(0,nY);
  
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
  
  # normalize the covariance by the variances yields the local family of correlation
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
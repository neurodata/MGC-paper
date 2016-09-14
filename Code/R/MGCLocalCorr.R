source("DistRanks.R")
source("DistCentering.R")

MGCLocalCorr <- function(X,Y,option, ind){
  # Author: Cencheng Shen
  # The main function that calculates all local correlation coefficients.
  #
  # The inputs are: 
  # two distance matrices X and Y;
  # an option that specifies which global correlation to use, including 'mcor','dcor','mantel'.
  #
  # The outputs are all local correlations and all local variances.
  #
  # Alternatively, specifying ind by a matrix single index will return a
  # weight matrix that shows the contribution of each distance entries to
  # the eventual local distance correlation at the given index.
  if (missing(option)){
    option='mcor'; # use mcorr by default
  }
  if (missing(ind)){
    ind=numeric(0);
  }
  n=nrow(X);
  disRank=cbind(DistRanks(X), DistRanks(Y)); # sort distances within columns
  
  # depending on the choice of the global correlation, properly center the distance matrices
  A=DistCentering(X,option);
  B=DistCentering(Y,option);
  
  RX=disRank[,1:n]; # the ranks for X
  RY=disRank[,(n+1):(2*n)]; # the ranks for Y
  if (length(ind)!=1){
    result=LocalCorrelations(A,t(B),RX,t(RY)); # compute all local corr / var statistics
  } else {
    result=LocalWeights(A,t(B),RX,t(RY),ind); # compute distance entry contributions at a given scale
  }
  return(result);
}

LocalCorrelations <- function(A,B,RX,RY){
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

LocalWeights <- function(A,B,RX,RY,ind){
  # An auxiliary function that computes the contributions of each distance entries to
  # the local distance correlation at a given scale.
  nX=max(RX);nY=max(RY);
  if (ind>nX*nY || ind<1){
    ind=nX*nY; # default to global scale when the specified index is out of range
  }
  k = ((ind-1) %% nX) + 1
  l = floor((ind-1) / nX) + 1
  RX=(RX>k);
  RY=(RY>l);
  A[RX]=0;
  B[RY]=0;
  weight=(A-mean(A))*(B-mean(B));
  return(weight);
}
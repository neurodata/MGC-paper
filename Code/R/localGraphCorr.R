source("disToRanks.R")

localGraphCorr <- function(X,Y,option, disRank){ # Calculate local graph correlation
  # Author: Cencheng Shen
  # Implements local graph correlation from Shen, Jovo, CEP 2016.
  # By specifying option=1, 2, or 3, it calculates the LGC statistics
  # by mcorr, dcorr, and Mantel
  if (missing(option)){
    option=1; # By default use mcorr
  }
  if (missing(disRank)){
    disRank=cbind(disToRanks(X), disToRanks(Y)); # Sort distances within columns, if the ranks are not provided
  }
  
  n=nrow(X);
  RX=disRank[,1:n]; # The ranks for X
  RY=disRank[,(n+1):(2*n)]; # The ranks for X
  corrXY=matrix(0,n,n);
  varX=rep(0,n);
  varY=rep(0,n);
  
  # Depending on the choice of global test, calculate the entries of A and B accordingly for late multiplication.
  if (option!=3){
    # Double centering for mcorr/dcorr
    H=diag(n)-(1/n)*matrix(1,n,n);
    A=H%*%X%*%H;
    B=H%*%Y%*%H;
    # For mcorr, further adjust the double centered matrices to remove high-dimensional bias
    if (option==1){
      A=A-X/n;
      B=B-Y/n;
      # The diagonals of mcorr are set to zero, instead of the original formulation of mcorr
      for (j in (1:n)){
        A[j,j]=0;
        B[j,j]=0;
      }
    }
  } else {
    # Single centering for Mantel
    EX=sum(X)/n/(n-1);
    EY=sum(Y)/n/(n-1);
    A=X-EX;
    B=Y-EY;
    # Mantel does not use diagonal entries, which is equivalent to set them zero
    for (j in (1:n)){
      A[j,j]=0;
      B[j,j]=0;
    }
  }
  
  # Summing up the entriwise product of A and B based on the ranks, which
  # yields the local family of covariance and variances
  for (j in (1:n)){
    for (i in (1:n)){
      a=A[i,j];
      b=B[i,j];
      # If there are ties, set all rank 0 entries to the diagonal entry
      if (RX[i,j]==0){
        a=A[j,j];
      }
      if (RY[i,j]==0){
        b=B[j,j];
      }
      tmp1=RX[i,j]+1;
      tmp2=RY[i,j]+1;
      corrXY[tmp1:n, tmp2:n]=corrXY[tmp1:n, tmp2:n]+a*b;
      varX[tmp1:n]=varX[tmp1:n]+a^2;
      varY[tmp2:n]=varY[tmp2:n]+b^2;
    }
  }
  
  # Normalizing the covariance by the variances yields the local family of correlation.
  options(warn=-1);
  corrXY=corrXY/Re(sqrt(varX%*%t(varY))); 
  options(warn=0);
  
  # Set any local correlation to 0 if any corresponding local variance is no larger than 0
  for (i in (1:n)){
    if (varX[i]<=0){
      corrXY[i,]=rep(0,n);
    }
    if (varY[i]<=0){
      corrXY[,i]=rep(0,n);
    }
  }
  result=list(corr=corrXY,varX=varX,varY=varY);
  return(result);
}
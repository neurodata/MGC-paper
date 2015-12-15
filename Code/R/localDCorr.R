source("disToRanks.R")

localDCorr <- function(X,Y,optionModified, disRank){ # Calculate local dCorr
  # Author: Cencheng Shen
  # Implements local distance correlation from Shen, Jovo, CEP 2015.
  if (missing(optionModified)){
    optionModified=1; # By default calculate the local version of both dcorr and mdcorr.
  }
  if (missing(disRank)){
    disRank=cbind(disToRanks(X), disToRanks(Y)); # Sort distances within columns, if the ranks are not provided
  }
  
  n=nrow(X);
  RX=disRank[,1:n];
  RY=disRank[,(n+1):(2*n)];
  corrXY=matrix(0,n,n);
  varX=rep(0,n);
  varY=rep(0,n);
  
  # Double centering of distance matrices
  H=diag(n)-(1/n)*matrix(1,n,n);
  XH=H%*%X%*%H;
  YH=H%*%Y%*%H;
  # Adjust the centering for modified dCorr
  if (optionModified==1){
    XH=n/(n-1)*(XH-X/n);
    YH=n/(n-1)*(YH-Y/n);
    for (i in (1:n)){
      XH[i,i]=0;
      YH[i,i]=0;
    }
  }
  
  # Summing up the product of distances, which yields dCov and dVar
  for (i in (1:n)){
    if (optionModified==0){
      corrXY=corrXY+XH[i,i]*YH[i,i];
      varX=varX+XH[i,i]*XH[i,i];
      varY=varY+YH[i,i]*YH[i,i];
    }
    for (j in (1:n)){
      if (i!=j){
        tmp1=RX[j,i]+1;
        tmp2=RY[j,i]+1;
        corrXY[tmp1:n, tmp2:n]=corrXY[tmp1:n, tmp2:n]+XH[j,i]*YH[j,i];
        varX[tmp1:n]=varX[tmp1:n]+XH[j,i]*XH[j,i];
        varY[tmp2:n]=varY[tmp2:n]+YH[j,i]*YH[j,i];
      }
    }
  }
  
  # Further adjust the diagonal products of modified dCov
  if (optionModified==1){
    meanX=sum(sum(X))/n^2;
    meanY=sum(sum(Y))/n^2;
    for (i in (1:n)){
      XH[i,i]=n/(n-1)*(mean(X[,i])-meanX);
      YH[i,i]=n/(n-1)*(mean(Y[,i])-meanY);
      corrXY=corrXY-2/(n-2)*XH[i,i]*YH[i,i];
      varX=varX-2/(n-2)*XH[i,i]*XH[i,i];
      varY=varY-2/(n-2)*YH[i,i]*YH[i,i];
    }
    corrXY=corrXY/n/(n-3);
    varX=varX/n/(n-3);
    varY=varY/n/(n-3);
  }
  # Normalizing dCov by dVar yields dCorr; there may exists NAN warning from below
  options(warn=-1);
  corrXY=corrXY/Re(sqrt(varX%*%t(varY))); 
  options(warn=0);
  
  # Set dCorr to 0 if any dVar is no larger than 0
  for (i in (1:(n-1))){
    if (varX[i]<=0){
      corrXY[i,]=rep(0,n);
    }
    if (varY[i]<=0){
      corrXY[,i]=rep(0,n);
    }
  }
  # The original dCorr is defined as the square root of the previous calculated dCorr; but square root or not does not affect testing at all.
  if (optionModified==0){
    corrXY=sqrt(corrXY);
  }
  result=list(corr=corrXY,varX=varX,varY=varY);
  return(result);
}
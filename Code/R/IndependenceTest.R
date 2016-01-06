source("disToRanks.R")
source("localGraphCorr.R")
library(ecodist)
library(HHG)
library(combinat)

IndependenceTest <-function(C,P,rep){
  # Author: Cencheng Shen
  # The independence Test for given data, an auxiliary function of the main
  # CorrPermDIstTest function
  if (missing(rep)){
    rep=1000; # By default use 1000 random permutations
  }
  
  n=nrow(C);
  alpha=0.05; # type 1 error level
  # Test statiscs under the null and the alternative
  dCor1N=array(0,dim=c(n,n,rep));dCor2N=array(0,dim=c(n,n,rep));dCor3N=array(0,dim=c(n,n,rep));
  dCor1A=array(0,dim=c(n,n,rep));dCor2A=array(0,dim=c(n,n,rep));dCor3A=array(0,dim=c(n,n,rep));
  # Powers for LGC by mcorr/dcorr/Mantel
  power1=matrix(0,n,n);power2=matrix(0,n,n);power3=matrix(0,n,n);
  
  # Calculate the rank matrices for given data 
  disRankC=disToRanks(C);
  disRankP=disToRanks(P);
  for (r in (1:rep)){
    # Random sampling with replacement
    per=sample(n,replace=TRUE);
    Ca=C[per,per];
    Pa=P[per,per];
    disRank=cbind(disRankC[per,per],disRankP[per, per]);
    # Calculate the test statistics under the alternative
    dCor1A[,,r]=localGraphCorr(Ca,Pa,1,disRank)$corr;
    dCor2A[,,r]=localGraphCorr(Ca,Pa,2,disRank)$corr;
    dCor3A[,,r]=localGraphCorr(Ca,Pa,3,disRank)$corr;
    
    # Random permutation
    perN=sample(n);
    perN=per[perN];
    Pa=P[perN,perN];
    disRank=cbind(disRankC[per,per],disRankP[perN, perN]);
    # Calculate the test statistics under the null
    dCor1N[,,r]=localGraphCorr(Ca,Pa,1,disRank)$corr;
    dCor2N[,,r]=localGraphCorr(Ca,Pa,2,disRank)$corr;
    dCor3N[,,r]=localGraphCorr(Ca,Pa,3,disRank)$corr;
  }
  
  # For each local test, estimate the critical value from the test statistics under the null,
  # then estimate the power from the test statistics under the alternative.
  for (i in (1:n)){
    for (j in (1:n)){
      dCorT=sort(dCor1N[i,j,],decreasing=TRUE);
      cut1=dCorT[ceiling(rep*alpha)];
      power1[i,j]=mean(dCor1A[i,j,]>cut1);
      
      dCorT=sort(dCor2N[i,j,],decreasing=TRUE);
      cut2=dCorT[ceiling(rep*alpha)];
      power2[i,j]=mean(dCor2A[i,j,]>cut2);
      
      dCorT=sort(dCor3N[i,j,],decreasing=TRUE);
      cut3=dCorT[ceiling(rep*alpha)];
      power3[i,j]=mean(dCor3A[i,j,]>cut3);
    }
  }
  
  # Set the powers of all local tests at rank 1 to 0 
  power1[1,]=rep(0,n);
  power2[1,]=rep(0,n);
  power3[1,]=rep(0,n);
  power1[,1]=rep(0,n);
  power2[,1]=rep(0,n);
  power3[,1]=rep(0,n);
  
  result=list(power1=power1,power2=power2,power3=power3);
  return(result);
}
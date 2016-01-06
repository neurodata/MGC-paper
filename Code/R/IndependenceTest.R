source("disToRanks.R")
source("localGraphCorr.R")
library(ecodist)
library(HHG)
library(combinat)

IndependenceTest <-function(C,P,rep){
  # Author: Cencheng Shen
  # Permutation Tests for identifying dependency, returning p-value of given data
  # The output are the p-values of local original dCorr, local modified dCorr, HHG, and Mantel test.
  
  # Parameters:
  # C and P are two distance matrices of same size for testing,
  # rep specifies the number of random permutations to use,
  # set allP to non-zero will use all permutations instead,
  # option specifies whether each test statistic is calculated or not.
  if (missing(rep)){
    rep=1000; # By default use 1000 random permutations
  }
  
  # P-values for local original dCorr, local modified dCorr, HHG, and Mantel test
  n=nrow(C);
  alpha=0.05;
  dCor1N=array(0,dim=c(n,n,rep));dCor2N=array(0,dim=c(n,n,rep));dCor3N=array(0,dim=c(n,n,rep));
  dCor1A=array(0,dim=c(n,n,rep));dCor2A=array(0,dim=c(n,n,rep));dCor3A=array(0,dim=c(n,n,rep));
  power1=matrix(0,n,n);power2=matrix(0,n,n);power3=matrix(0,n,n);
  
  # Calculate the test statistics for the given data sets
  disRankC=disToRanks(C);
  disRankP=disToRanks(P);
  # Now Permute the second dataset for rep times, and calculate the p-values
  for (r in (1:rep)){
    # Use random permutations; if allP is not 0, use all possible permutations
    per=sample(n,replace=TRUE);
    Ca=C[per,per];
    Pa=P[per,per];
    disRank=cbind(disRankC[per,per],disRankP[per, per]);
    dCor1A[,,r]=localGraphCorr(Ca,Pa,1,disRank)$corr;
    dCor2A[,,r]=localGraphCorr(Ca,Pa,2,disRank)$corr;
    dCor3A[,,r]=localGraphCorr(Ca,Pa,3,disRank)$corr;
    
    perN=sample(n);
    Pa=P[perN,perN];
    disRank=cbind(disRankC[per,per],disRankP[perN, perN]);
    dCor1N[,,r]=localGraphCorr(Ca,Pa,1,disRank)$corr;
    dCor2N[,,r]=localGraphCorr(Ca,Pa,2,disRank)$corr;
    dCor3N[,,r]=localGraphCorr(Ca,Pa,3,disRank)$corr;
  }
  
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
  
  # Treat the p-value of local methods in neighborhood 0 as 1
  power1[1,]=rep(0,n);
  power2[1,]=rep(0,n);
  power3[1,]=rep(0,n);
  power1[,1]=rep(0,n);
  power2[,1]=rep(0,n);
  power3[,1]=rep(0,n);
  
  result=list(power1=power1,power2=power2,power3=power3);
  return(result);
}
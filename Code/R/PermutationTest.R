source("disToRanks.R")
source("localGraphCorr.R")
library(ecodist)
library(HHG)
library(combinat)

PermutationTest <-function(C,P,rep,allP,option){
  # Author: Cencheng Shen
  # The permutation Test for given data, an auxiliary function of the main
  # CorrPermDIstTest function
  if (missing(rep)){
    rep=1000;
  }
  if (missing(allP)){
    allP=0; 
  }
  if (missing(option)){
    option=c(1,1,1,1); 
  }
  if (allP!=0){
    PALL=permn(n); 
    rep=length(PALL);
  }
  
  # P-values for LGC by mcorr, LGC by dcorr, LGC by Mantel, and HHG
  n=nrow(C);
  p1=matrix(0,n,n); p2=matrix(0,n,n);p3=matrix(0,n,n);p4=0;
  
  # Calculate the observed test statistics for the given data sets
  disRankC=disToRanks(C);
  disRankP=disToRanks(P);
  disRank=cbind(disRankC,disRankP);
  if (option[1]!=0){
    cut1=localGraphCorr(C,P,1,disRank)$corr;
  }
  if (option[2]!=0){
    cut2=localGraphCorr(C,P,2,disRank)$corr;
  }
  if (option[3]!=0){
    cut3=localGraphCorr(C,P,3,disRank)$corr;
  }
  if (option[4]!=0){
    cut4=hhg.test(C,P,nr.perm=0);
    cut4=unlist(cut4[1]);
  }
  
  # Now Permute the second dataset for rep times, and calculate the p-values
  for (r2 in (1:rep)){
    # Use random permutations; if allP is not 0, use all possible permutations
    per=sample(n);
    if (allP!=0){
      per=unlist(PAll[r2]);
    }
    Pa=P[per,per];
    disRank=cbind(disRankC,disRankP[per, per]);
    if (option[1]!=0){
      dCor1=localGraphCorr(C,Pa,1,disRank)$corr;
      p1=p1+(dCor1<cut1)*1/rep;
    }
    if (option[2]!=0){
      dCor2=localGraphCorr(C,Pa,2,disRank)$corr;
      p2=p2+(dCor2<cut2)*1/rep;
    }
    if (option[3]!=0){
      dCor3=localGraphCorr(C,Pa,3,disRank)$corr;
      p3=p3+(dCor3<cut3)*1/rep;
    }
    if (option[4]!=0){
      dCor4=hhg.test(C,Pa,nr.perm=0);
      dCor4=unlist(dCor4[1]);
      p4=p4+(dCor4<cut4)*1/rep;
    }
  }
  
  # Output the p-value
  p1=1-p1;
  p2=1-p2;
  p3=1-p3;
  p4=1-p4;
  
  result=list(LGCmcorr=p1,LGCdcorr=p2,LGCMantel=p3,HHG=p4,mcorr=p1[n,n],dcorr=p2[n,n],Mantel=p3[n,n]);
  return(result);
}
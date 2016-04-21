source("disToRanks.R")
source("LocalGraphCorr.R")
#library(ecodist)
#library(combinat)

MGCPermutationTest <-function(C,D,rep,option){
  # Author: Cencheng Shen
  # This is the main function for MGC permutation test, based on
  # dcorr/mcorr/Mantel.
  # The output are estimated MGC p-value, all p-values of local test, the
  # estimated optimal scale, the estimated MGC test statistic, and all local
  # test statistics.
  
  # Parameters:
  # C & D should be two n*n distance matrices,
  # rep specifies the number of random permutations to use for the permutation test,
  # option specifies which global test to use.
  if (missing(rep)){
    rep=1000;
  }
  if (missing(option)){
    option=1; # Default option. Set to 1/2/3 for local dcorr/mcorr/Mantel
  }
  n=nrow(C);
  pAll=matrix(0,n,n); 
  
  # Calculate the observed test statistics for the given data sets
  testAll=LocalGraphCorr(C,D,option)$corr;
  
  # Now Permute the second dataset for rep times, and calculate the p-values
  for (r in (1:rep)){
    # Use random permutations; 
    per=sample(n);
    DN=D[per,per];
    tmp=LocalGraphCorr(C,DN,1)$corr;
    pAll=pAll+(tmp<testAll)*1/rep;
  }
  
  # Set the p-values of rank 0 to maximum.
  pAll=1-pAll;
  pAll[1,1:n]=1;pAll[1:n,1]=1;
  
  # Verify and estimate the MGC optimal scale
  ind=Verify(pAll);
  # Get the estimated MGC p-value and MGC test statistic
  if (length(ind)==0){
    p=1;
    test=0;
  }else{
    p=pAll[ind];
    test=testAll[ind];
  }
  
  result=list(pMGC=p,pAll=pAll,MGC=test,MGCAll=testAll,ind=ind);
  return(result);
}

Verify <-function(pAll){
  n=nrow(pAll);
  alpha=0.05;thres=0.85;delta=n/20;
  power1=(pAll<alpha)*1;
  pCol=colMeans(power1[2:n,1:n]);
  pRow=rowMeans(power1[1:n,2:n]);
  k=which(pCol>thres);
  l=which(pRow>thres);
  col=FALSE;
  
  if (length(k)<length(l)){
    pAll=t(pAll);
    pRow=pCol;
    k=l;
    col=TRUE;
  }
  if (length(k)<delta){
    ind=numeric();
    if (pRow[n]>thres || length(which(pRow>pRow[n]))==0){
      ind=n^2;
    }
    return(ind);
  }
  
  k=which(pRow==max(pRow));
  ind=which(pAll[k,]==min(pAll[k,]));
  ind=ind[length(ind)];
  
  l = ceiling(ind / length(k)); 
  k1 = ind -(l-1)*length(k);
  k=k[k1];
  if (col==TRUE){
    tmp=l;
    l=k;
    k=tmp;
  }
  ind=(k-1)*n+l;
  print(k)
  print(l)
  print(ind)
}
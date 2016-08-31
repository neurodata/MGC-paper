source("disToRanks.R")
source("LocalCorr.R")
source("MGCScaleVerify.R")

MGCPermutationTest <-function(A,B,rep,option){
  # Author: Cencheng Shen
  # This function tests independent between two data sets, using MGC by a random permutation test.
  #
  # The inputs are:
  # two n*n distance matrices A & B,
  # a parameter rep to specify the number of random permutations,
  # an option to specify which global test to use, set to 1,2,3 for dcorr / mcorr / Mantel.
  #
  # The outputs are:
  # the estimated MGC p-value, the p-values of all local tests, 
  # the estimated MGC test statistic, all local test statistics, and the estimated optimal scale. 
  if (missing(rep)){
    rep=1000; # use 1000 random permutations by default
  }
  if (missing(option)){
    option=1; # use mcorr by default
  }
  n=nrow(A);
  
  # calculate all local correlations between the two data sets
  testAll=LocalCorr(A,B,option)$corr;
  
  # calculate the local correlations under permutation, to yield the p-values of all observed local correlations
  for (r in (1:rep)){
    # use random permutations on the second data set
    per=sample(n);
    BN=B[per,per];
    tmp=LocalCorr(A,BN,option)$corr;
    if (r==1){
      pAll=(tmp<testAll)*1/rep;
    }else{
      pAll=pAll+(tmp<testAll)*1/rep;
    }
  }
  
  # set the p-values of local corr at rank 0 to maximum, since they should not be used
  pAll=1-pAll;
  pAll[1,]=1;pAll[,1]=1;
  
  # verify and estimate the MGC optimal scale
  result=MGCScaleVerify(pAll);
  p=result$p;
  ind=result$ind;
  test=testAll[ind[length(ind)]];
  
  result=list(pMGC=p,pAll=pAll,MGC=test,testAll=testAll,ind=ind);
  return(result);
}
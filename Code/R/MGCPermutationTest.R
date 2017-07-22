source("MGCLocalCorr.R")
source("MGCSampleStat.R")
source("FindLargestRectangles.R")

MGCPermutationTest <-function(A,B,rep,option){
  # Author: Cencheng Shen
  # This function tests independent between two data sets, using MGC by a random permutation test.
  # It outputs all local correlations, the sample MGC statistic, their
  # p-values, and the estimated optimal scales.
  #
  # The inputs are:
  # two distance matrices A & B,
  # a parameter rep to specify the number of random permutations,
  # an option to specify which global test to use, by 'mcor', 'dcor', 'mantel'.
  #
  # The outputs are:
  # the sample MGC p-value, sample MGC test statistic,
  # all local p-values, all local correlations, and the estimated optimal scales.
  #
  # Note that the sample MGC p-value / test statistic / optimal scales are
  # not calculated unless the global correlation is 'specified as mcor', and
  # the optimal scales are output as matrix single indices.
  if (missing(rep)){
    rep=1000; # use 1000 random permutations by default
  }
  if (missing(option)){
    option='mcor'; # use mcorr by default
  }
  
  sampleIndicator=0;
  if (option=='mcor'||option=='mgc'){
    sampleIndicator=1; # only compute sample MGC for mcorr
  } 
  # calculate all local correlations between the two data sets
  localCorr=MGCLocalCorr(A,B,option)$corr;
  if (sampleIndicator==1){
    statMGC=MGCSampleStat(localCorr) # sample MGC for the observed data
  } 
  m=nrow(localCorr);
  n=ncol(localCorr);
  pLocalCorr=matrix(0,m,n);
  pMGC=0;
  n2=nrow(B);
  
  # calculate the local correlations under permutation, to yield the p-values of all observed local correlations
  for (r in (1:rep)){
    # use random permutations on the second data set
    per=sample(n2);
    BN=B[per,per];
    tmp=MGCLocalCorr(A,BN,option)$corr;
    pLocalCorr=pLocalCorr+(tmp>=localCorr)*1/rep;
    if (sampleIndicator==1){
      tmp2=MGCSampleStat(tmp) # sample MGC for permuted data
      pMGC=pMGC+(tmp2>=statMGC)*1/rep;
    }
  }
  if (sampleIndicator!=1){ # other than mcorr, we do not implemented sample MGC yet, and the global statistic is always used
    pMGC=pLocalCorr[m,n];
    statMGC=localCorr[m,n];
  }
  
  # if p-value equals 0, enlarge it to 1/rep, since the identity permutation is always
  # one such permutation.
  # if (min(pLocalCorr)==0){
    # pLocalCorr=pLocalCorr+1/rep;
  # }
  # pLocalCorr[pLocalCorr>1]=1;
  # pLocalCorr[1,]=1;pLocalCorr[,1]=1;
  # if (min(pLocalCorr[2:m,2:n])>pMGC){
    # pMGC=min(pLocalCorr[2:m,2:n]);
  # }
  
  # estimate the optimal scales
  optimalInd=FindLargestRectangles((pLocalCorr<=pMGC)&(localCorr>=statMGC))$M;
  optimalInd=which(optimalInd==1);
  if (pLocalCorr[m,n]<pMGC || length(which(optimalInd==m*n))>0){
    optimalInd=m*n; # if the global scale is not selected in the largest rectangle while being optimal, we take the global scale instead.
  }
  
  result=list(pMGC=pMGC,statMGC=statMGC,pLocalCorr=pLocalCorr,localCorr=localCorr,optimalInd=optimalInd);
  return(result);
}
source("disToRanks.R")
source("LocalCorr.R")

MGCPermutationTest <-function(A,B,rep,option){
  # Author: Cencheng Shen
  # This function implements a random permutation test for MGC.
  # The output are estimated MGC p-value, the p-values of all local tests, the estimated MGC test statistic,
  # all local test statistics, and the estimated optimal scale, 
  
  # Parameters:
  # A & B should be two n*n distance matrices,
  # rep specifies the number of random permutations to use for the permutation test,
  # option specifies which global test to use.
  if (missing(rep)){
    rep=1000;
  }
  if (missing(option)){
    option=1; # default option, set to 1/2/3 for local dcorr/mcorr/Mantel
  }
  n=nrow(A);
  
  # calculate the observed test statistics for the given data sets
  testAll=LocalCorr(A,B,option)$corr;
  
  # permute the second dataset for rep times, and calculate the p-values
  for (r in (1:rep)){
    # use random permutations; 
    per=sample(n);
    BN=B[per,per];
    tmp=LocalCorr(A,BN,option)$corr;
    if (r==1){
      pAll=(tmp<testAll)*1/rep;
    }else{
      pAll=pAll+(tmp<testAll)*1/rep;
    }
  }
  
  # set the p-values of rank 0 to maximum.
  pAll=1-pAll;
  pAll[1,]=1;pAll[,1]=1;
  
  # verify and estimate the MGC optimal scale
  ind=MGCScaleVerify(pAll);
  p=pAll[ind];
  test=testAll[ind];
  
  result=list(pMGC=p,pAll=pAll,MGC=test,testAll=testAll,ind=ind);
  return(result);
}

MGCScaleVerify <-function(V){
  # An auxiliary function to verify and estimate the MGC optimal scale
  VN=V[2:nrow(V),2:ncol(V)];
  k=Verify(VN)+1;
  l=Verify(t(VN))+1;
  ind=which(V==min(V[k,l]));
  return(ind[length(ind)]);
}

Verify <-function(VN){
  thres=0.05;
  m=nrow(VN);
  k=m;
  rowTmp=rep(0,m);
  rowTmp2=rep(0,m);
  for (i in (1:m)){
    rowTmp2[i]=median(VN[i,]);
    rowTmp[i]=min(VN[i,]);
  }
  indK=which(rowTmp2==min(rowTmp2));
  indK=indK[length(indK)];
  if (rowTmp2[indK]<=thres){
    rowTmp=(rowTmp<thres);
    tmp=indK;
    for (i in ((indK-1):1)){
      if (i < 1 || rowTmp[i]==FALSE){
        break;
      } else {
        tmp=append(tmp,i);
      }
    }
    for (i in ((indK+1):m)){
      if (i > m || rowTmp[i]==FALSE){
        break;
      } else {
        tmp=append(tmp,i);
      }
    }
    VN=VN[tmp,];
    if (median(VN)<=thres/m*length(tmp)){
      k=tmp;
    }
  }
  return(k);
}
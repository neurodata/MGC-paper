library(ecodist)
library(SDMTools)
source("MGCLocalCorr.R")

MGCSampleStat <- function(A,B,option){
  # Author: Cencheng Shen
  # This function estimate the Oracle MGC (i.e., optimal local correlation)
  # from the local correlation map, which we call sample MGC statistic.
  #
  # It finds the largest connected region in the correlation map, such that
  # each correlation is significant, i.e., larger than a certain threshold.
  # To avoid correlation inflation by sample noise, it then computes Sample MGC as follows:
  # for the largest correlation in the region, calcualte the two minimal correlations
  # along adjacent row scales and adjacent column scales, then take the larger one as the Sample MGC.
  # If the region area is too small, or the estimated Sample MGC is no larger than the global, use the global correlation instead.
  #
  # Input: either a size m*n local correlation map, or two n*n distance matrices A and B and the global corr option.
  # Output: the sample MGC statistic within [-1,1].
  if (missing(B) && missing(option)){
    localCorr=A; # if there is only one input, asume the localCorr is given as A
  } else {
    localCorr=MGCLocalCorr(A,B,option)$corr; # otherwise compute the localCorr from given distance matrices
  }
  m=nrow(localCorr);
  n=ncol(localCorr);

  negCorr=localCorr[2:m,2:n];
  negCorr=negCorr[negCorr<0]; # negative correlations
  thres1=sqrt(sum(negCorr^2)/length(negCorr));  # threshold based on negative correlations
  if (is.na(thres1) || thres1<0.01){
    thres1=0.01;
  }
  thres1=thres1*3.5;
  thres2=min(2/min(m,n),0.05); # threshold based on sample size

  statMGC=localCorr[m,n]; # take the global correlation by default

  R=(localCorr>max(thres1,thres2)); # find all correlations that are larger than the threshold
  # find the largest connected component of all significant correlations
  if (sum(R)>0){
     R=ConnCompLabel(R==1);
     tmp=tabulate(R);
     tmp=which.max(tmp);
     R=(R==tmp);
  }

  if (mean(R)>=thres2){ # proceed only when the region area is sufficiently large
      ind=which((localCorr>=max(localCorr[R]))&(R==1)); # find the scale within R that has the maximum correlation
      k = ((ind-1) %% m) + 1
      l = floor((ind-1) / m) + 1

      ln=ceiling(0.1*m); # boundary for checking adjacent rows
      km=ceiling(0.1*n); # boundary for checking adjacent columns
      for (i in (1:length(k))){
        ki=k[i];
        li=l[i];
        
        # ensure the adjacent rows does not exceed the local correlation size, same for columns
        left=max(2,li-ln); 
        right=min(n,li+ln);
        upper=max(2,ki-km);
        down=min(m,ki+km);
        
        tmp1=min(localCorr[upper:down,li]); # minimal correlation at given row and adjacent columns
        tmp2=min(localCorr[ki,left:right]); # minimal correlation at given column and adjacent rows
        tmp=max(tmp1,tmp2); # take the max of the two minimal correlations
        if (tmp>statMGC){
          statMGC=tmp; 
        }
      }
  }
  return(statMGC);
}
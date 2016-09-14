library(ecodist)
library(SDMTools)
MGCSampleStat <- function(localCorr){
  # Author: Cencheng Shen
  # This function estimate the Oracle MGC (i.e., optimal local correlation)
  # from the local correlation map, which we call sample MGC statistic.
  #
  # It finds a sufficiently large region in the correlation map, such that
  # all correlations within the region are significantly large and
  # monotonically changing. Then a relatively large correlation in the region
  # is used as sample MGC. If no such region exists, use the global
  # correlation instead.
  #
  # Input: a size m*n local correlation map.
  # Output: the sample MGC statistic within [-1,1].
  m=nrow(localCorr);
  n=ncol(localCorr);
  negCorr=localCorr[2:m,2:n];
  negCorr=negCorr[negCorr<0]; # negative correlations
  
  eps=sqrt(sum(negCorr^2)/length(negCorr));
  if (is.na(eps) || eps<0.01){
    eps=0.01;
  }
  eps=3.5*eps; # threshold for significantly large correlation
  
  R=(localCorr>eps); # find a region with significantly large correlation
  thres=min(2/min(m,n),0.05); # threshold for sufficiently large region
  statMGC=localCorr[m,n]; # take global correlation by default

  if (mean(R)>=thres){
    R=Monotone(localCorr,R); # put monotonically changing restriction to reduce the region R
  
    if (mean(R)>0){
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
        tmp=max(tmp1,tmp2); # take the max of the minimal correlations, which should be relatively large but does not inflate much over the true optimal statistic
        if (tmp>statMGC){
          statMGC=tmp; 
        }
      }
    }
  }
  return(statMGC);
}

Monotone <- function(localCorr,R){
  # An auxiliary function that computes regions in R with monotonically increasing
  # or decreasing correlations along the row or column
  m=nrow(localCorr);
  n=ncol(localCorr);
  
  PD1=matrix(0,m,n);
  PD2=matrix(0,m,n);
  for (i in (2:m)){
    tt=localCorr[i,];
    PD1[i,2:n]=diff(tt); # store the p-value changes within rows
  }
  for (i in (2:n)){
    tt=localCorr[,i];
    PD2[2:m,i]=diff(tt); # store the p-value changes within columns
  }
  
  R2=list(matrix(0,m,n),matrix(0,m,n),matrix(0,m,n),matrix(0,m,n));
  R2[[1]]=(PD1>=0)&R; # monotonically increasing along the row within R
  R2[[2]]=(PD1<=0)&R; # monotonically decreasing along the row within R
  R2[[3]]=(PD2>=0)&R; # monotonically increasing along the column within R
  R2[[4]]=(PD2<=0)&R; # monotonically decreasing along the column within R
  R=as.logical(matrix(0,m,n));
  
  for (i in (1:4)){
    t=sum(R2[[i]]);
    if (t>0){
      t=ConnCompLabel(R2[[i]]); # connected component of each region 
      R= (R | t); # combine all monotonically changing and significant regions
    }
    return(R)
  }
}
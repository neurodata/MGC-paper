source("FindLargestRectangles.R")

MGCScaleVerify <-function(P,gamma,tau){
  # Author: Cencheng Shen
  # This function approximates the p-value and the optimal scale for sample MGC,
  # based on the p-value map of all local correlations.
  # The algorithm first searches for large smooth regions in the p-value map.
  # If the global p-value is small enough among all local p-values, the global p-value is used for MGC;
  # or if there exists a smooth rectangular region that is large enough,
  # a smaller p-value within the smooth rectangle is approximated for MGC.
  # Once we determine the sample MGC p-value, the smooth rectangle that is
  # bounded above by the p-value is taken as the optimal scales, which shall include the global scale if necessary.
  # Otherwise, we use the mean p-value of all local p-values instead.
  if (missing(gamma)){
    gamma=0.1; # gamma is used to: determine if the global p-value is significant enough, determine if the rectangular region is significant enough, and approximate a small p-value in the significant rectangular region
  }
  if (missing(tau)){
    tau=0.005; # tau is a threshold to approximate the monotone p-values change for smooth regions
  }
  m=nrow(P);
  n=ncol(P);
  gamma=max(gamma,3/min(m,n)); # increase gamma accordingly in case the sample size is too small
  
  # find large smooth regions in the p-value map
  R=SmoothRegions(P,tau);
  
  # check for global p-value and smooth rectangle
  if (mean(P<P[m,n])<=gamma){
    p=P[m,n]; # directly use the global p-value if it is among the top 100*gamma% of all local p-values
  } else {
    R=FindLargestRectangles(R)$M; # find the largest rectangle within the smooth regions
    tmp=mean(R);
    # approximate a smaller p-value from the smooth rectangle if and only if the area is no smaller than gamma
    if (tmp>=gamma){
      # take a small p-value that is 100*gamma/2area(R)% of all p-values in the smooth rectangle. 
      # For example, if the area of R equals gamma, the median
      # p-value within R is use for the MGC p-value.
      p=quantile(P[R],gamma/2/tmp);
    }else{
      p=mean(P[P<1]); # otherwise, use the mean p-value for MGC
    }
  }
  
  # the largest smooth rectangle bounded by the sample MGC p-value is taken as the optimal scale
  R=FindLargestRectangles((R==1) & (P<=p))$M;
  if (sum(R)==0){
    # if the scales from above return empty, relax the smooth region contraint for optimal scale. 
    # This rarely happens except when the mean p-value is taken.
    R=FindLargestRectangles((P<=p))$M;
  }
  
  # include the global scale if necessary
  ind=which(R==1);
  if (p>=P[m*n] && R[m*n]==0){
    ind=c(ind,m*n);
  }
  
  result=list(p=p,ind=ind);
  return(result);
}

SmoothRegions <- function(P,tau){
  # Find the smooth regions in the p-value map, by considering the
  # largest monotonically decreasing or increasing scales along
  # the row or column p-values, but allowing small p-value increase or
  # decrease as specified by tau. 
  m=nrow(P);
  n=ncol(P);
  
  PD=list(matrix(0,m,n),matrix(0,m,n));
  for (i in (2:m)){
    tt=P[i,];
    PD[[1]][i,2:n]=diff(tt); # store the p-value changes within rows
  }
  for (i in (2:n)){
    tt=P[,i];
    PD[[2]][2:m,i]=diff(tt); # store the p-value changes within columns
  }
  
  RC=list(matrix(0,m,n),matrix(0,m,n),matrix(0,m,n),matrix(0,m,n));
  RC[[1]]=(PD[[1]]<=tau); # check monotone decreasing , but also allows small p-value increase no more than tau
  RC[[2]]=(PD[[1]]>=-tau); # check monotone increasing , but also allows small p-value decrease no more than tau
  RC[[3]]=(PD[[2]]<=tau); # repeat for column changes
  RC[[4]]=(PD[[2]]>=-tau);
  
  for (i in (1:4)){
    RC[[i]]=FindLargestRectangles(RC[[i]])$M; # find the largest rectangle within each (approximately) monotonically decreasing / increasing region
  }
  
  R=(RC[[1]] | RC[[2]] | RC[[3]] | RC[[4]]); # combine all four rectangles together into the smooth regions
}
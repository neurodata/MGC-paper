MGCScaleVerify <-function(V){
  # An auxiliary function to verify and estimate the MGC optimal scale based
  # on the p-values of all local correlations.
  VN=V[2:nrow(V),2:ncol(V)];
  k=Verify(VN)+1;# verify the rows
  l=Verify(t(VN))+1;# verify the columns
  ind=which(V==min(V[k,l]));
  return(ind[length(ind)]);
}

Verify <-function(VN){
  thres=0.05; # type 1 error level
  m=nrow(VN);
  k=m; # take the largest scale by default
  
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
    
    # from the row with minimal median p-value, include adjacency rows
    # whose median p-value are significant
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
    
    # check if the included rows satisfy a Bonferroni-type bound
    VN=VN[tmp,];
    if (median(VN)<=thres/m*length(tmp)){
      k=tmp;
    }
  }
  return(k);
}
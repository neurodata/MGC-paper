library(ggplot2)
library(reshape2)
library(lattice)
source("PermutationTest.R")
source("IndependenceTest.R")
source("verifyNeighbors.R")

CorrPermDistTest <- function(type,rep,cv,titlechar,allP,option){
  # Author: Cencheng Shen
  # Permutation Tests for identifying dependency.
  # The output are the p-values of LGC by mcorr/dcorr/Mantel, and HHG;
  # followed by the estimated optimal neighborhood for LGC by
  # mcorr/dcorr/Mantel.
  # Note that the local family include the global test at the last entry.
  
  # Parameters:
  # type should be a n*2n matrix, a concatenation of two distance matrices,
  # rep specifies the number of random permutations to use,
  # cv specifies the number of bootstrap samples to use for neighborhood validation,
  # set allP to non-zero will use all permutations instead,
  # option specifies whether each test statistic is calculated or not.
  if (missing(rep)){
    rep=1000; # Default random permutation numbers
  }
  if (missing(allP)){
    allP=0; # If set to other value, will override rep and use all permutations; unfeasible for n large
  }
  if (missing(titlechar)){
    titlechar="Data";
  }
  if (missing(cv)){
    cv=1000; # Default bootstrap replicates to estimate the optimal neighborhood 
  }
  if (missing(option)){
    option=c(1,1,1,1);  # Default option. Setting any to 0 to disable the calculation of mcorr/dcorr/Mantel/HHG.
  }
  
  n=nrow(type);
  C=type[, 1:n];
  P=type[, (n+1):(2*n)];
  
  # If cv is not 0, use resampling to estimate the optimal neighborhood by the testing powers
  ps1=matrix(0,n,n);ps2=matrix(0,n,n);
  if (cv!=0){
    testP=IndependenceTest(C,P,cv);
    neighbor1=verifyNeighbors(1-testP$power1);
    neighbor2=verifyNeighbors(1-testP$power2);
    neighbor3=verifyNeighbors(1-testP$power3);
  }
  # Return p-values from the permutation test
  testP=PermutationTest(C,P,rep,allP,option);
  if (cv==0){
    neighbor1=verifyNeighbors(testP$p1);
    neighbor2=verifyNeighbors(testP$p2);
    neighbor3=verifyNeighbors(testP$p3);
  }
  
  # Plot level plot
  max=0.2;p=testP$LGCmcorr;
  p[which(p>max)]=max;
  myAt=seq(0,max,0.02);
  interval=5;
  ckey=list(at=myAt,labels=list(cex=2));
  col.l <- colorRampPalette(c('red', 'orange', 'yellow', 'green', 'cyan', 'blue'))
  ckey=list(labels=list(cex=2));
  lplot=levelplot(p,zscaleLog="e",col.regions = terrain.colors(100),at=myAt,scales=list(x=list(at=seq(interval,n,interval), cex=2), y=list(at=seq(interval,n,interval), cex=2)),xlab=list(label="Neighborhood Choice of X",cex=2),ylab=list(label="Neighborhood Choice of Y",cex=2),main=list(label="Permutation Test P-Value",cex=2),colorkey=ckey)
  
  # Output
  output=list(titlechar=titlechar,LGCmcorr=mean(testP$LGCmcorr[neighbor1]),LGCdcorr=mean(testP$LGCdcorr[neighbor2]),LGCMantel=mean(testP$LGCMantel[neighbor3]),HHG=testP$HHG, mcorr=testP$mcorr,dcorr=testP$dcorr,Mantel=testP$Mantel,n=n,rep=rep,allP=allP,option=option,levelPlot=lplot);
  #output=list(titlechar=titlechar,LGCmcorr=testP$LGCmcorr,LGCdcorr=testP$LGCdcorr,LGCMantel=testP$LGCMantel,HHG=testP$HHG, mcorr=testP$mcorr,dcorr=testP$dcorr,Mantel=testP$Mantel,n=n,rep=rep,allP=allP,option=option,neighbor1=neighbor1,neighbor2=neighbor2,neighbor3=neighbor3);
  
  filename = paste("CorrPermDistTestType", titlechar,".RData",sep = "");
  save(output,file = filename);
  return(output);
}
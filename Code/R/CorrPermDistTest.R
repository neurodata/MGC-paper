library(ggplot2)
library(reshape2)
source("PermutationTest.R")
source("verifyNeighbors.R")

CorrPermDistTest <- function(type,rep,cv,titlechar,allP,option){
  # Author: Cencheng Shen
  # Permutation Tests for identifying dependency, returning p-value of given data
  # The output are the p-values of local original dCorr, local modified dCorr, HHG, and Mantel test.
  
  # Parameters:
  # type should be a n*2n matrix, a concatenation of two distance matrices,
  # rep specifies the number of random permutations to use,
  # cv specifies the number of bootstrap samples to use for neighborhood validation
  # set allP to non-zero will use all permutations instead,
  # option specifies whether each test statistic is calculated or not.
  if (missing(rep)){
    rep=1000; # By default use 1000 random permutations
  }
  if (missing(allP)){
    allP=0; # By default do not use all permutations
  }
  if (missing(titlechar)){
    titlechar="Data";
  }
  if (missing(cv)){
    cv=0;
  }
  if (missing(option)){
    option=c(1,1,1,1); # Control whether to calculate the respective correlation statistic or not.
  }
  
  n=nrow(type);
  C=type[, 1:n];
  P=type[, (n+1):(2*n)];
  
  ps1=matrix(0,n,n);ps2=matrix(0,n,n);
  if (cv!=0){
    ratio=0.5;
    optionCV=option;
    optionCV[3:4]=rep(0,2);
    for (i in (1:cv)){
      noise=rnorm(n,0,1);
      noise=as.matrix(dist(noise));
      noise=noise/norm(noise,'f')*norm(P,'f')*ratio;
      Pa=P+noise;
      #per=sample(n,n,replace=TRUE);
      per=(1:n);
      testP=PermutationTest(C[per,per],Pa[per,per],rep,allP,optionCV);
      ps1=ps1+testP$ldcorr/cv;
      ps2=ps2+testP$lmdcorr/cv;
    }
  }
  testP=PermutationTest(C,P,rep,allP,option);
  if (cv==0){
    ps1=ps1+testP$ldcorr;
    ps2=ps2+testP$lmdcorr;
  }
  neighbor1=verifyNeighbors(ps1);
  neighbor2=verifyNeighbors(ps2);
  
  output=list(titlechar=titlechar,ldcorr=testP$ldcorr,lmdcorr=testP$lmdcorr,HHG=testP$HHG,Mantel=testP$Mantel,dcorr=testP$dcorr,mdcorr=testP$mdcorr,dcorrNeighbor=neighbor1,mdcorrNeighbor=neighbor2,n=n,rep=rep,allP=allP,option=option);
  
  # Plot the p-value w.r.t. neighborhood
  n=output$n;
  p1=rep(0,n);p2=rep(0,n);p3=rep(output$dcorr,n);p4=rep(output$mdcorr,n);p5=rep(output$HHG[1],n);p6=rep(output$Mantel,n);
  for (i in (1:n)){
     p1[i]=min(output$ldcorr[i,]);
     p2[i]=min(output$lmdcorr[i,]);
  }
  plotData=data.frame(LDcorr=p1,LMDcorr=p2,DCorr=p3,MDCorr=p4,HHG=p5,Mantel=p6,Neighborhood=(1:n));
  plotData <- melt(plotData,id="Neighborhood");
  colnames(plotData) <- c("Neighborhood","Method","value");
  p=ggplot(data=plotData,aes(x=Neighborhood,y=value,color=Method,linetype=Method))+geom_line()+geom_point()+scale_y_log10()+scale_linetype_manual(values = c(rep("solid", 2), rep("dashed", 4)))+scale_colour_manual(values=c("blue","red","blue","red","green","cyan"))+ylab("P-Value")+ggtitle("Permutation Test for Brain Data")+theme(plot.title = element_text(size=20,face="bold"),legend.title = element_text(size=16, face="bold"),legend.text = element_text(size=12),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"));

  
  # Save the results
  filename = paste("CorrPermDistTestType", titlechar,".RData",sep = "");
  save(output,file = filename);
  return(output);
}
### dcorr Code validation by Brain CxP, to check out dcorr is the same as the energy pack.
rm(list = ls())
load("BrainCPData.RData")
load("BrainHippoShape.RData")
source("localGraphcorr.R")
library(ecodist)
library(energy)
library(HHG)
C=distC;P=distP;
C=LMRS;P=as.matrix(dist(Label))+1;
for (i in (1:n)){
  P[i,i]=0;
}
# Our own local dcorr and local mdcorr
ldcorr=localDCorr(C,P,0)$corr;
lmdcorr=localDCorr(C,P,1)$corr;
# Global dcorr from energy package.
# The output squared should be the same as last entry in ldcorr
dcorr=dcor(as.dist(C),as.dist(P));
# Global mdcorr from energy package.
# The bias correted dcor in the output is the same as the last entry in lmdcorr. 
mdcorr=dcor.ttest(C,P,distance=TRUE);
# Mantel test
mantel=mantel(as.dist(C)~as.dist(P),nperm=1000);
# HHG test
hhgr=hhg.test(C,P,nr.perm=1000);


### dCorr permutation test on the brain data
rm(list = ls())
load("BrainCPData.RData")
C=distC;P=distP;
source("CorrPermDistTest.R")
test=CorrPermDistTest(cbind(C,P),rep=100,cv=100)
min(min(test$ldcorr))
min(min(test$lmdcorr))
test$HHG
test$Mantel

load("BrainHippoShape.RData")
C=LMLS;P=as.matrix(dist(Label));
test=CorrPermDistTest(cbind(C,P),rep=1000)
min(min(test$ldcorr))
min(min(test$lmdcorr))
test$HHG
test$Mantel

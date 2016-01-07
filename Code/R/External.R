### dCorr permutation test on the brain data
rm(list = ls())
load("BrainCPData.RData")
C=distC;P=distP;
source("CorrPermDistTest.R")
test=CorrPermDistTest(cbind(C,P),rep=10000,cv=10000,"BrainCxP")
LGCmcorr=mean(testP$LGCmcorr[testP$neighbormcorr]);
LGCdcorr=mean(testP$LGCdcorr[testP$neighbordcorr]);
LGCMantel=mean(testP$LGCMantel[testP$neighborMantel]);

rm(list = ls())
load("BrainHippoShape.RData")
C=LMLS;P=as.matrix(dist(Label))+1;diag(P)=0;
test=CorrPermDistTest(cbind(C,P),rep=1000,cv=1000,"BrainLMLY")
C=LMRS;
test=CorrPermDistTest(cbind(C,P),rep=1000,cv=1000,"BrainLMRY")

### dcorr Code validation by Brain CxP, to check out dcorr is the same as the energy pack.
rm(list = ls())
load("BrainCPData.RData")
source("localGraphcorr.R")
library(ecodist)
library(energy)
library(HHG)
C=distC;P=distP;
# Our own local dcorr and local mdcorr
ldcorr=localGraphCorr(C,P,0)$corr;
lmdcorr=localGraphCorr(C,P,1)$corr;
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

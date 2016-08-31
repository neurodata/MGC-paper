# load the data
load("BrainHippoShape.RData")
### Compute all local corr
library(ecodist)
library(energy)
library(HHG)
source("LocalCorr.R")
C=LMLS;P=as.matrix(dist(Label))+1;diag(P)=0;
ldcorr=LocalCorr(C,P,option=1)$corr;
lmdcorr=LocalCorr(C,P,option=2)$corr;
lmantel=LocalCorr(C,P,option=3)$corr;
### Global corr from the original authors. 
### Their p-values in permutation test should be very similar to our local corr at largest scale,
### but the actual statistic can be slightly different.
# global dcorr
dcorr=dcor(as.dist(C),as.dist(P));
# global mcorr
mdcorr=dcor.ttest(C,P,distance=TRUE);
# Mantel test
mantel=mantel(as.dist(C)~as.dist(P),nperm=1000);
# HHG test
hhgr=hhg.test(C,P,nr.perm=1000);


### Permutation Test of local corr
source("MGCPermutationTest.R")
test=MGCPermutationTest(C,P,rep=100,option=2)
C=LMRS;
test=MGCPermutationTest(C,P,rep=100,option=2)



load("BrainCPData.RData")
source("LocalCorr.R")
library(ecodist)
library(energy)
library(HHG)
C=distC;P=distP;
ldcorr=LocalCorr(C,P,option=1)$corr;
lmdcorr=LocalCorr(C,P,option=2)$corr;
lmantel=LocalCorr(C,P,option=3)$corr;
dcorr=dcor(as.dist(C),as.dist(P));
mdcorr=dcor.ttest(C,P,distance=TRUE);
mantel=mantel(as.dist(C)~as.dist(P),nperm=1000);
hhgr=hhg.test(C,P,nr.perm=1000);

### Permutation Test of local corr
source("MGCPermutationTest.R")
test=MGCPermutationTest(C,P,rep=1000,option=2)

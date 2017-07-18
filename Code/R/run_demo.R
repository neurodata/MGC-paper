library(ecodist)
library(energy) #dcorr package
library(HHG) #hhg package

load("../../Data/Preprocessed/BrainCPData.RData")


C=distC;P=distP;
mdcorr=dcor.ttest(C,P,distance=TRUE); #unbiased dcorr test
mantel=mantel(as.dist(C)~as.dist(P),nperm=1000); #mantel test
hhgr=hhg.test(C,P,nr.perm=1000); # hhg test

source("MGCPermutationTest.R")
mgc=MGCPermutationTest(C,P,rep=1000,option='mcor') #mgc test

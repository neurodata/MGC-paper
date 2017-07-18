%% simulate data
type=6;
d=1;
noise=0;
n=100;

[x,y]=CorrSampleGenerator(type,n,d,0, noise);


%% run mgc
rep=1000;
option='mcor';

A=squareform(pdist(x));
B=squareform(pdist(y));

[pMGC,statMGC,pLocalCorr,localCorr,optimalInd]=MGCPermutationTest(A,B,rep,option);

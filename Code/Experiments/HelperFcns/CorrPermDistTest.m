function [pMGC,pD,pM,pP, pHHG]=CorrPermDistTest(C,D,rep,titlechar,option)
% Author: Cencheng Shen
% Permutation Tests for identifying dependency.
% The output are the p-values of MGC by dcorr/mcorr/Mantel, and global dcorr/mcorr/Mantel/HHG.

% Parameters:
% C & D should be two n*n distance matrices,
% rep specifies the number of random permutations to use for the permutation test,
% alpha specifies the type 1 error level,
% option specifies whether each test statistic is calculated or not,
% neighborhood can be specified beforehand so as to skip the optimal scale
% estimation.
if nargin<4
    titlechar='RealData';
end
if nargin<5
    option=[1,2,3,4]; % Default option. Setting any to 0 to disable the calculation of MGC by dcorr/mcorr/Mantel, global dcorr/mcorr/Mantel, HHG, in order; set the first three
end

% Global Correlations
pDLocal=1;pMLocal=1;pPLocal=1;pD=1;pM=1;pP=1;optimalInd=0;testMGC=0;pMGC=1;testDLocal=0;testMLocal=0;testPLocal=0;pHHG=1;testHHG=0;
if option(1)==1
    [~,~,pDLocal,testDLocal,~]=MGCPermutationTest(C,D,rep,'dcor');
end
if option(2)==2
    [pMGC,testMGC,pMLocal,testMLocal,optimalInd]=MGCPermutationTest(C,D,rep,'mcor');
end
if option(3)==3
    [~,~,pPLocal,testPLocal,~]=MGCPermutationTest(C,D,rep,'mantel');
end
if option(4)==4
    [pHHG,testHHG]=HHGPermutationTest(C,D,rep);
end
pD=pDLocal(end);pM=pMLocal(end);pP=pPLocal(end);

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-3));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

pre1=strcat(rootDir,'Data/Results/');% The folder to save figures
filename=strcat(pre1,'CorrPermDistTestType',titlechar);
save(filename,'titlechar','rep','option','pDLocal','pMLocal','pPLocal','pHHG','pMGC','pD','pM','pP','optimalInd','testMGC','testDLocal','testMLocal','testMLocal','testHHG');

function  [pHHG, testHHG]=HHGPermutationTest(C,D,rep)
% Author: Cencheng Shen
% This is for HHG permutation test.
% The outputs are the p-values of HHG and the test statistic.
if nargin<3
    rep=1000;
end
n=size(C,1);

% Calculate the observed test statistics for the given data sets
pHHG=0;
testHHG=HHG(C,D);

% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % Use random permutations;
    per=randperm(n);
    DN=D(per,per);
    tmp=HHG(C,DN);
    pHHG=pHHG+(tmp>=testHHG)/rep;
end
if (sum(pHHG<1/rep)>0)
    pHHG=pHHG+1/rep;
end
function  [pMGC,testMGC,pLocal,testLocal,optimalInd]=MGCPermutationTest(A,B,rep,option)
% Author: Cencheng Shen
% This function tests independent between two data sets, using MGC by a random permutation test.
%
% The inputs are:
% two n*n distance matrices A & B,
% a parameter rep to specify the number of random permutations,
% an option to specify which global test to use, set to 1,2,3 for mcor / dcor / Mantel.
%
% The outputs are:
% the estimated MGC p-value, the p-values of all local tests,
% the estimated MGC test statistic, all local test statistics, and the estimated optimal scale.

if nargin<3
    rep=1000; % use 1000 random permutations by default
end
if nargin<4
    option='mcor';  % use mcorr by default
end
[m,n]=size(A);

if strcmp(option,'mcor')
    ind=1;
else
    ind=0;
end
% calculate all local correlations between the two data sets
testLocal=LocalCorr(A,B,option);
if ind==1
    testMGC=SampleMGC(testLocal);
end
pLocal=zeros(size(testLocal));pMGC=0;
% calculate the local correlations under permutation, to yield the p-values of all observed local correlations
for r=1:rep
    % use random permutations on the second data set
    per=randperm(n);
    BN=B(per,per);
    tmp=LocalCorr(A,BN,option);
    tmp2=SampleMGC(tmp);
    pLocal=pLocal+(tmp>=testLocal)/rep;
    if ind==1
        pMGC=pMGC+(tmp2>=testMGC)/rep;
    end
end
if ind~=1
    pMGC=1;
    testMGC=0;
end
% set the p-values of local corr at rank 0 to maximum, since they should not be used
if min(min(pLocal(2:end,2:end)))==0
    pLocal=pLocal+1/rep;
end
if min(min(pLocal(2:end,2:end)))>pMGC
    pMGC=min(min(pLocal(2:end,2:end)));
end
pLocal(pLocal>1)=1;
pLocal(1,:)=1;pLocal(:,1)=1;

% find optimal scale
warning('off','all');
[~,~,~,optimalInd]=FindLargestRectangles((pLocal<=pMGC), [0 0 1],[2,2]);
optimalInd=find(optimalInd==1);
if (pLocal(end)<pMGC && isempty(find(optimalInd==m*n, 1)))
    optimalInd=m*n;
end
function  [pMGC,testMGC,pLocal,testLocal,optimalInd]=MGCPermutationTest(A,B,rep,option,alpha)
% Author: Cencheng Shen
% This function tests independent between two data sets, using MGC by a random permutation test.
%
% The inputs are:
% two n*n distance matrices A & B,
% a parameter rep to specify the number of random permutations,
% an option to specify which global test to use, set to 1,2,3 for dcorr / mcorr / Mantel.
%
% The outputs are:
% the estimated MGC p-value, the p-values of all local tests,
% the estimated MGC test statistic, all local test statistics, and the estimated optimal scale.

if nargin<3
    rep=1000; % use 1000 random permutations by default
end
if nargin<4
    option=2;  % use mcorr by default
end
if nargin<5
    alpha=0.05;  % use mcorr by default
end
[m,n]=size(A);

% calculate all local correlations between the two data sets
testLocal=LocalCorr(A,B,option);
if option==2
    testMGC=SampleMGC(testLocal);
end
pLocal=0;pMGC=0;
% calculate the local correlations under permutation, to yield the p-values of all observed local correlations
for r=1:rep
    % use random permutations on the second data set
    per=randperm(n);
    BN=B(per,per);
    tmp=LocalCorr(A,BN,option);
    tmp2=SampleMGC(tmp);
    pLocal=pLocal+(tmp>=testLocal)/rep;
    if option==2
        pMGC=pMGC+(tmp2>=testMGC)/rep;
    end
end
% set the p-values of local corr at rank 0 to maximum, since they should not be used
if pMGC==0 || min(min(pLocal(2:end,2:end)))==0
    pLocal=pLocal+1/rep;
    pMGC=pMGC+1/rep;
end
pLocal(pLocal>1)=1;
pLocal(1,:)=1;pLocal(:,1)=1;

% find optimal scale
if option==2
    warning('off','all');
    if min(min(pLocal))>pMGC
        pMGC=min(min(pLocal));
    end
    [~,~,~,optimalInd]=FindLargestRectangles((pLocal<=pMGC), [0 0 1],[2,2]);
    optimalInd=find(optimalInd==1);
    if pLocal(end)<=pMGC && (isempty(optimalInd) || sum(optimalInd==m*n)==0)
        optimalInd=[optimalInd;m*n];
    end
else
    testMGC=testLocal(end);
    pMGC=pLocal(end);
    optimalInd=m*n;
end

if pMGC>alpha
    optimalInd=1;
end
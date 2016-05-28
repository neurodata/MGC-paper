function  [p,pAll,test,testAll,ind]=MGCPermutationTest(A,B,rep,option)
% Author: Cencheng Shen
% This function implements a random permutation test for MGC.
% The output are estimated MGC p-value, the p-values of all local tests, the estimated MGC test statistic,
% all local test statistics, and the estimated optimal scale, 

% Parameters:
% A & B should be two n*n distance matrices,
% rep specifies the number of random permutations to use for the permutation test,
% option specifies which global test to use.
if nargin<3
    rep=1000;
end
if nargin<4
    option=1;  % default option, set to 1/2/3 for local dcorr/mcorr/Mantel
end
n=size(A,1);

% calculate the observed test statistics for the given data sets
testAll=LocalCorr(A,B,option);

% permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % use random permutations;
    per=randperm(n);
    BN=B(per,per);
    tmp=LocalCorr(A,BN,option);
    if r==1
        pAll=(tmp<testAll)/rep;
    else
        pAll=pAll+(tmp<testAll)/rep;
    end
end
% set the p-values of rank 0 to maximum.
pAll=1-pAll;
pAll(1,:)=1;pAll(:,1)=1;

% verify and estimate the MGC optimal scale
ind=MGCScaleVerify(pAll);
p=pAll(ind);
test=testAll(ind);
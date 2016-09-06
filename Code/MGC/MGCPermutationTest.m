function  [p,pAll,test,testAll,indAll]=MGCPermutationTest(A,B,rep,option)
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
[m,n]=size(A);

% calculate all local correlations between the two data sets
testAll=LocalCorr(A,B,option);
 test=SampleMGC(testAll);
% calculate the local correlations under permutation, to yield the p-values of all observed local correlations
for r=1:rep
    % use random permutations on the second data set
    per=randperm(n);
    BN=B(per,per);
    tmp=LocalCorr(A,BN,option);
    tmp2=SampleMGC(tmp);
    if r==1
        pAll=(tmp>=testAll)/rep;
        p=(tmp2>=test)/rep;
    else
        pAll=pAll+(tmp>=testAll)/rep;
        p=p+(tmp2>=test)/rep;
    end
end
% set the p-values of local corr at rank 0 to maximum, since they should not be used
if (sum(sum(pAll<1/rep))>0)
    pAll=pAll+1/rep;
end
pAll(pAll>1)=1;
pAll(1,:)=1;pAll(:,1)=1;

% find optimal scale
warning('off','all');
[~,~,~,indAll]=FindLargestRectangles((pAll<=p), [0 0 1],[2,2]);
indAll=find(indAll==1);
[m,n]=size(pAll);
if pAll(end)<=p && sum(indAll==m*n)==0
    indAll=[indAll;m*n];
end
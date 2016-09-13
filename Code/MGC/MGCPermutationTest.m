function  [pMGC,statMGC,pLocalCorr,localCorr,optimalInd]=MGCPermutationTest(A,B,rep,option)
% Author: Cencheng Shen
% This function tests independent between two data sets, using MGC by a random permutation test.
% It outputs all local correlations, the sample MGC statistic, their
% p-values, and the estimated optimal scales.
%
% The inputs are:
% two distance matrices A & B,
% a parameter rep to specify the number of random permutations,
% an option to specify which global test to use, by 'mcor', 'dcor', 'mantel'.
%
% The outputs are:
% the sample MGC p-value, sample MGC test statistic,
% all local p-values, all local correlations, and the estimated optimal scales.
%
% Note that the sample MGC p-value / test statistic / optimal scales are
% not calculated unless the global correlation is 'specified as mcor', and
% the optimal scales are output as matrix single indices.

if nargin<3
    rep=1000; % use 1000 random permutations by default
end
if nargin<4
    option='mcor';  % use mcorr by default
end

if strcmp(option,'mcor')
    ind=1; % only compute sample MGC for mcorr
else
    ind=0;
end

% calculate all local correlations between the two data sets
localCorr=LocalCorr(A,B,option);
[m,n]=size(localCorr);
if ind==1
    statMGC=MGCSampleStat(localCorr); % sample MGC for the observed data
end
pLocalCorr=zeros(size(localCorr));pMGC=0;
n2=size(B,2);

% calculate the local correlations under permutation, to yield the p-values of all observed local correlations
for r=1:rep
    % use random permutations on the second data set
    per=randperm(n2);
    BN=B(per,per);
    tmp=LocalCorr(A,BN,option);
    tmp2=MGCSampleStat(tmp);
    pLocalCorr=pLocalCorr+(tmp>=localCorr)/rep;
    if ind==1
        pMGC=pMGC+(tmp2>=statMGC)/rep;
    end
end
if ind~=1
    pMGC=1;
    statMGC=0;
end
% if p-value equals 0, enlarge it to 1/rep, since the identity permutation is always
% one such permutation.
if min(min(pLocalCorr(2:end,2:end)))==0
    pLocalCorr=pLocalCorr+1/rep;
end
if min(min(pLocalCorr(2:end,2:end)))>pMGC
    pMGC=min(min(pLocalCorr(2:end,2:end)));
end
pLocalCorr(pLocalCorr>1)=1;
pLocalCorr(1,:)=1;pLocalCorr(:,1)=1;

% estimate the optimal scales
warning('off','all');
[~,~,~,optimalInd]=FindLargestRectangles((pLocalCorr<=pMGC), [0 0 1],[2,2]);
optimalInd=find(optimalInd==1);
if (pLocalCorr(end)<pMGC && isempty(find(optimalInd==m*n, 1)))
    optimalInd=m*n; % if the global scale is not selected in the largest rectangle while being optimal, we take the global scale instead.
end
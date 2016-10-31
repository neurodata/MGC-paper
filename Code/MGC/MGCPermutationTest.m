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

sampleIndicator=0;
if strcmp(option,'mcor')
    sampleIndicator=1; % only compute sample MGC for mcorr
end

% calculate all local correlations between the two data sets
localCorr=MGCLocalCorr(A,B,option);
[m,n]=size(localCorr);
n2=size(B,2);

tmp=zeros(100,m,n);
for r=1:100
    % use random permutations on the second data set
    n2=size(B,2);
    per=randperm(n2);
    BN=B(per,per);
    tmp(r,:,:)=MGCLocalCorr(A,BN,option);
%     if sampleIndicator==1
%         tmp(r)=MGCSampleStat(tmp2,0.03); % sample MGC for permuted data
%     end
end
thres=zeros(m,n);
for i=2:m
    for j=2:n
        thres(i,j)=prctile(tmp(:,i,j),95);
    end
end
    [k,l]=find(thres>=max(max(thres))); % find the scale within R that has the maximum correlation
    
    thresM=0;
    ln=ceil(0.1*m); % boundary for checking adjacent rows
    km=ceil(0.1*n); % boundary for checking adjacent columns
    for i=1:length(k)
        ki=k(i);
        li=l(i);
        
        % ensure the adjacent rows does not exceed the local correlation size, same for columns
        left=max(2,li-ln);
        right=min(n,li+ln);
        upper=max(2,ki-km);
        down=min(m,ki+km);
        
        tmp1=min(thres(upper:down,li)); % minimal correlation at given row and adjacent columns
        tmp2=min(thres(ki,left:right)); % minimal correlation at given column and adjacent rows
        tmp=max(tmp1,tmp2); % take the max of the two minimal correlations
        if tmp>thresM
            thresM=tmp;
        end
    end
thres=thresM


if sampleIndicator==1
    statMGC=MGCSampleStat(localCorr,thres); % sample MGC for the observed data
end
pLocalCorr=zeros(size(localCorr));pMGC=0;


% calculate the local correlations under permutation, to yield the p-values of all observed local correlations
for r=1:rep
    % use random permutations on the second data set
    per=randperm(n2);
    BN=B(per,per);
    tmp=MGCLocalCorr(A,BN,option);
    pLocalCorr=pLocalCorr+(tmp>=localCorr)/rep;
    if sampleIndicator==1
        tmp2=MGCSampleStat(tmp,thres); % sample MGC for permuted data
        pMGC=pMGC+(tmp2>=statMGC)/rep;
    end
end
if sampleIndicator~=1 % other than mcorr, we do not implemented sample MGC yet, and the global statistic is always used
    pMGC=pLocalCorr(end);
    statMGC=localCorr(end);
end
% if p-value equals 0, enlarge it to 1/rep, since the identity permutation is always
% one such permutation.
if min(min(pLocalCorr(2:end,2:end)))==0
    pLocalCorr=pLocalCorr+1/rep;
end
pLocalCorr(pLocalCorr>1)=1;
pLocalCorr(1,:)=1;pLocalCorr(:,1)=1;
if min(min(pLocalCorr(2:end,2:end)))>pMGC
    pMGC=min(min(pLocalCorr(2:end,2:end)));
end

% estimate the optimal scales
[~,~,~,optimalInd]=FindLargestRectangles((pLocalCorr<=pMGC), [0 0 1],[2,2]);
optimalInd=find(optimalInd==1);
if (pLocalCorr(end)<=pMGC && isempty(find(optimalInd==m*n, 1)))
    optimalInd=m*n; % if the global scale is not selected in the largest rectangle while being optimal, we take the global scale instead.
end
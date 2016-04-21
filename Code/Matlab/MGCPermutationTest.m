function  [p,pAll,test,testAll,ind]=MGCPermutationTest(C,D,rep,option)
% Author: Cencheng Shen
% This is the main function for MGC permutation test, based on
% dcorr/mcorr/Mantel.
% The output are estimated MGC p-value, all p-values of local test, the
% estimated optimal scale, the estimated MGC test statistic, and all local
% test statistics.

% Parameters:
% C & D should be two n*n distance matrices,
% rep specifies the number of random permutations to use for the permutation test,
% option specifies which global test to use.
if nargin<3
    rep=1000;
end
if nargin<4
    option=1;  % Default option. Set to 1/2/3 for local dcorr/mcorr/Mantel
end
n=size(C,1);
pAll=zeros(n,n);

% Calculate the observed test statistics for the given data sets
testAll=LocalGraphCorr(C,D,option);

% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % Use random permutations;
    per=randperm(n);
    DN=D(per,per);
    tmp=LocalGraphCorr(C,DN,option);
    pAll=pAll+(tmp<testAll)/rep;
end
% Set the p-values of rank 0 to maximum.
pAll=1-pAll;
pAll(1,:)=1;pAll(:,1)=1;

% Verify and estimate the MGC optimal scale
ind=Verify(pAll);

% Output the estimated MGC p-value and MGC test statistic
if isempty(ind)
    p=1;
    test=0;
else
    p=pAll(ind);
    test=testAll(ind);
end

function ind=Verify(pAll)
% An auxiliary function to verify and estimate the MGC optimal scale
n=size(pAll,1);
alpha=0.05;
thres=0.85;
delta=n/20;
power1=(pAll<alpha);
pCol=mean(power1(2:end,1:end),1);
pRow=mean(power1(1:end,2:end),2);
k=find(pRow>thres);
l=find(pCol>thres);
col=false;

if length(k)<length(l)
    pAll=pAll';
    pRow=pCol;
    k=l;
    col=true;
end
if length(k)<delta;% || mean(a(indK))>length(indK)
    ind=[];
    if pRow(n)>thres || sum(pRow>pRow(n))==0
        ind=n^2;
    end
    return;
end
% figure
% plot(1:n,pCol,'b-',1:n,pRow,'r:')

k=find(pRow==max(pRow));
ind=find(pAll(k,:)==min(min(pAll(k,:))),1,'last');
[k1,l]=ind2sub([length(k),n],ind);
k=k(k1);
if col==true
    tmp=l;
    l=k;
    k=tmp;
end
ind=sub2ind([n,n],k,l);


function [p,ind]=verify3(pAll)
% alpha=0.06;
thres=0.1;
if pAll==1;
    p=1;
    ind=[];
    return;
end
n=size(pAll,1);
% power1=(p1All<alpha);
pCol=mean(pAll(2:end,1:end),1);
pRow=mean(pAll(1:end,2:end),2);

% figure
% plot(1:n,pCol,'b-',1:n,pRow,'r:')

if min(pRow)>min(pCol)
    pAll=pAll';
    pRow=pCol;
    pCol=pRow;
end
% [~,~,a]=unique(1-pRow);
% k=find(a==1,1,'last');
k=find(pRow<thres);
if length(k)<n/10;% || mean(a(indK))>length(indK)
    ind=[];
    p=1;
    if pRow(n)<thres || sum(pRow<pRow(n))==0
        ind=n^2;
        p=pAll(ind);
    end
    return;
end
k=k(end);
l=find(pAll(k,:)==min(pAll(k,:)),1,'last');
p=pAll(k,l);
if min(pRow)>min(pCol)
    tmp=l;
    l=k;
    k=tmp;
end
ind=sub2ind([n,n],k,l);
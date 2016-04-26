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
% pAll=zeros(n,n);

% Calculate the observed test statistics for the given data sets
testAll=LocalGraphCorr(C,D,option);

% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % Use random permutations;
    per=randperm(n);
    DN=D(per,per);
    tmp=LocalGraphCorr(C,DN,option);
    if r==1
        pAll=(tmp<testAll)/rep;
    else
        pAll=pAll+(tmp<testAll)/rep;
    end
end
% Set the p-values of rank 0 to maximum.
pAll=1-pAll;
pAll(1,:)=1;pAll(:,1)=1;

% Verify and estimate the MGC optimal scale
ind=MGCScaleVerify(pAll);
p=pAll(ind);
test=testAll(ind);
    
% % Output the estimated MGC p-value and MGC test statistic
% if isempty(ind)
%     p=1;
%     test=0;
% else
%     p=pAll(ind);
%     test=testAll(ind);
% end

% function [ind]=Verify(pAll)
% % An auxiliary function to verify and estimate the MGC optimal scale
% [nX,nY]=size(pAll);
% thres=0.5;
% delta=0;%1/50;
% alpha=0.05;
% gamma=0.02;
% % power1=(pAll<alpha);
% % pCol=mean(power1(2:end,1:end),1);
% % pRow=mean(power1(1:end,2:end),2);
% pRow=zeros(nX,1);pCol=zeros(nY,1);
% for i=1:nX
%     pRow(i)=mean(pAll(i,2:end)<min(min(pAll(i,2:end)),alpha)+gamma);  
% end
% for i=1:nY
%     pCol(i)=mean(pAll(2:end,i)<min(min(pAll(2:end,i)),alpha)+gamma);  
% end
% 
% 
% k=find(pRow>thres);
% l=find(pCol>thres);
% rowB=nX;
% if length(k)<length(l)
%     pAll=pAll';
%     pRow=pCol;
%     k=l;
%     rowB=nY;
% end
% 
% if isempty(k);
%     ind=nX*nY;
%     return;
% end
% % figure
% % plot(2:nX,pRow(2:end));
% 
% k=find(pRow==max(pRow));
% ind=find(pAll(k,:)==min(min(pAll(k,:))),1,'last');
% [k1,l]=ind2sub(size(pAll(k,:)),ind);
% k=k(k1);
% % pAll(k,l)
% if rowB==nY
%     tmp=l;
%     l=k;
%     k=tmp;
%     pAll=pAll';
% end
% ind=sub2ind(size(pAll),k,l);
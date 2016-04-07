function [p1,p2,p3, p4,p5,p6,p7,p1All,p2All,p3All]=CorrPermDistTest(C,D,rep1,rep2,titlechar,option)
% Author: Cencheng Shen
% Permutation Tests for identifying dependency.
% The output are the p-values of MGC by dcorr/mcorr/Mantel, and global dcorr/mcorr/Mantel/HHG.

% Parameters:
% C & D should be two n*n distance matrices,
% rep1 specifies the number of independent samples for optimal scale estimation
% in MGC; if 0, the estimation step is skipped.
% rep2 specifies the number of random permutations to use for the permutation test,
% alpha specifies the type 1 error level,
% option specifies whether each test statistic is calculated or not,
% neighborhood can be specified beforehand so as to skip the optimal scale
% estimation.
if nargin<5
    titlechar='RealData';
end
if nargin<6
    option=[1,2,3,4]; % Default option. Setting any to 0 to disable the calculation of MGC by dcorr/mcorr/Mantel, global dcorr/mcorr/Mantel, HHG, in order; set the first three
end
n=size(C,1);
ratio1=0.5;
nn=floor(n*ratio1);

% Global Correlations
[p1All]=PermutationTest(C,D,rep2,option(1));
[p2All]=PermutationTest(C,D,rep2,option(2));
[p3All]=PermutationTest(C,D,rep2,option(3));
[p7]=PermutationTest(C,D,rep2,option(4));
p4=p1All(end);
p5=p2All(end);
p6=p3All(end);
if rep1==0;
    ind1=maxNeighbors(1-p1All);ind2=maxNeighbors(1-p2All);ind3=maxNeighbors(1-p3All);
else
    % % Run the permutation test to return the p-values
    if rep2>1000
        rr=floor(rep2/rep1);
    else
        rr=rep2;
    end
    ind1=splitPermutation(C,D,nn,rep1,rr,option(1));
    ind2=splitPermutation(C,D,nn,rep1,rr,option(2));
    ind3=splitPermutation(C,D,nn,rep1,rr,option(3));
end
% p1=median(p1All(ind1));
% p2=median(p2All(ind2));
% p3=median(p3All(ind3));
p1=min(p1All(ind1));
p2=min(p2All(ind2));
p3=min(p3All(ind3));
if min(p1,p4)<0.05
    p1
    p4
end
% Save the results
pre1='../../Data/';
filename=strcat(pre1,'CorrPermDistTestType',titlechar);
save(filename,'titlechar','rep2','option','p1All','p2All','p3All','p4','p5','p6','p7','p1','p2','p3','ind1','ind2','ind3');

function ind=splitPermutation(C,D,nn,rep1,rep2,option)
ind=[];
if option==0;
    return;
end
n=size(C,1);
epsilon=0.2;
ind=1:nn^2;
for r=1:rep1
    per=randperm(n);
    l1=length(ind);
    for rr=1:floor(n/nn)
        pa1=per((rr-1)*nn+1:rr*nn);
        [tmp]=PermutationTest(C(pa1,pa1),D(pa1,pa1),rep2,option);
        ind=ind(tmp(ind)-min(min(tmp))<epsilon);
    end
    length(ind)
    if length(ind)<5;
        ind=[];
        break;
    end
    if length(ind)==l1
        break;
    end
end
if isempty(ind);
    ind=n^2;
else
    ind=ratioScale(ind,nn,n);
end

function ind=splitPermutation2(C,D,nn,rep1,rep2,option)
ind=[];
if option==0;
    return;
end
n=size(C,1);
for r=1:rep1
    per=randperm(n);
    pMean=zeros(nn,nn);
    %     pStd=ones(nn,nn);
    %dCorA=ones(nn,nn,floor(n/nn));
    for rr=1:floor(n/nn)
        pa1=per((rr-1)*nn+1:rr*nn);
        [tmp]=PermutationTest(C(pa1,pa1),D(pa1,pa1),rep2,option);
        %         tmp=tmp-min(min(tmp(2:end,2:end)));
        pMean=pMean+tmp/floor(n/nn);
    end
    %     for i=2:nn
    %         for j=2:nn
    %             tmp=reshape(dCorA(i,j,:),1,floor(n/nn));
    %             pMean(i,j)=mean(tmp);
    % %             pStd(i,j)=std(tmp);
    %         end
    %     end
    
    %     crStd=sort(pStd(2:end,2:end),'ascend');
    %     crStd=crStd(ceil(length(crStd)*0.3));
    if r==1
        ind=find(pMean-min(min(pMean))<0.05);
    else
        ind=ind(pMean(ind)-min(min(pMean))<0.1);
    end
    length(ind)
    %     ind=ind(pStd(ind)<0.2);
    if isempty(ind)
        length(ind)
        break
    end
end
% length(ind)
if isempty(ind);
    length(ind)
    ind=n^2;
else
    ind=ratioScale(ind,nn,n);
end

function [ind]=ratioScale(indn,nn,n)
% nn=size(p1All,1);
% indn=maxNeighbors(1-p1All,0);
% indn=maxNeighbors(p1All,0);
indLength=length(indn);
ind=[];
% ratio=(n-1)/(nn-1);
for i=1:indLength;
    [k,l]=ind2sub([nn,nn],indn(i));
    k=k-1;
    l=l-1;
    if k==0 || l==0
        continue;
    end
    kb=floor((k-1)*(n-1)/(nn-1))+2;
    kt=floor(k*(n-1)/(nn-1))+1;
    lb=floor((l-1)*(n-1)/(nn-1))+2;
    lt=floor(l*(n-1)/(nn-1))+1;
    for kk=kb:kt
        for ll=lb:lt
            tmp=sub2ind([n,n], kk, ll);
            ind=[ind tmp];
        end
    end
end

function  [p1]=PermutationTest(C,D,rep,option)
% This is an auxiliary function of the main function to calculate the p-values of
% all local tests of dcorr/mcorr/Mantel, the p-value of HHG in the
% permutation test.
if nargin<4
    option=1;  % Default option. Set to 1/2/3/4 for local dcorr/mcorr/Mantel, or HHG.
end

n=size(C,1);
if option==0
    p1=1;
    return;
end

% Calculate the observed test statistics for the given data sets
if option~=4
    p1=zeros(n,n);
    cut1=LocalGraphCorr(C,D,option);
    dCor1=zeros(n,n,rep);  % Default option. Set to 1/2/3/4 for local dcorr/mcorr/Mantel, or HHG.
else
    p1=0;
    cut1=HHG(C,D);
    dCor1=0;
end

% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % Use random permutations;
    per=randperm(n);
    DN=D(per,per);
    CN=C;
    if option<4
        dCorTmp=LocalGraphCorr(CN,DN,option);
        dCor1(:,:,r)=dCorTmp;
    else
        dCorTmp=HHG(CN,DN);
        dCor1(r)=dCorTmp;
    end
    p1=p1+(dCorTmp<cut1)/rep;
end
% Output the p-value
p1=1-p1;
if option~=4
    p1(1,:)=1;p1(:,1)=1; % Set the p-values of rank 0 to maximum.
end
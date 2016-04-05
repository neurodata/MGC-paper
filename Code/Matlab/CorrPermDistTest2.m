function [p1,p2,p3, p4,p5,p6,p7,p1All,p2All,p3All]=CorrPermDistTest2(C,D,rep1,rep2,titlechar,option)
% Author: Cencheng Shen
% Permutation Tests for identifying dependency.
% The output are the p-values of MGC by mcorr/dcorr/Mantel, and global mcorr/dcorr/Mantel/HHG.

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
    option=[1,2,0,0]; % Default option. Setting any to 0 to disable the calculation of MGC by mcorr/dcorr/Mantel, global mcorr/dcorr/Mantel, HHG, in order; set the first three
end
ratio=0.5;
n=size(C,1);
nn=ceil(n*ratio);
% if nn<10
%     nn=10;
% end
p1=0;p2=0;p3=0;

if rep1~=0
    % % Run the permutation test to return the p-values
    p1=splitPermutation(C,D,nn,rep1,rep2,option(1));
    p2=splitPermutation(C,D,nn,rep1,rep2,option(2));
    p3=splitPermutation(C,D,nn,rep1,rep2,option(3));
end

% Global Correlations
[p1All]=PermutationTest(C,D,rep2,option(1));
[p2All]=PermutationTest(C,D,rep2,option(2));
[p3All]=PermutationTest(C,D,rep2,option(3));
[p7]=PermutationTest(C,D,rep2,option(4));
p4=p1All(end);p5=p2All(end);p6=p3All(end);
if rep1==0;
    %ind1=maxNeighbors(1-p1All);ind2=maxNeighbors(1-p2All);ind3=maxNeighbors(1-p3All);
    p1=min(min(p1All));
    p2=min(min(p2All));
    p3=min(min(p3All));
end
% Save the results
pre1='../../Data/';
filename=strcat(pre1,'CorrPermDistTestType',titlechar);
save(filename,'titlechar','rep2','option','p1All','p2All','p3All','p4','p5','p6','p7','p1','p2','p3');

function p1=splitPermutation(C,D,nn,rep1,rep2,option)
if nargin<5
    option=1;
end
if option==0;
    p1=1;
    return;
end
p1=zeros(rep1,1);
dCor1=0;
n=size(C,1);
for r=1:rep1
    per=randperm(n);
    pa1=per(1:nn);
    pa2=per(nn+1:n);
    [p1All]=PermutationTest(C(pa1,pa1),D(pa1,pa1),rep2,option);
    ind1=ratioScale(p1All,n);
    [p1All]=PermutationTest(C(pa2,pa2),D(pa2,pa2),rep2,option);
    p1(r)=median(p1All(ind1));
end
%     if rep1>=3
%         rem=ceil(rep1*0.1);
%         p1=sort(p1);p1=p1(rem+1:rep1-rem);
%     end
p1=mean(p1);

function ind=ratioScale(p1All,n)
nn=size(p1All,1);
p1All(1,:)=1;p1All(:,1)=1;
indn=maxNeighbors(1-p1All);
indLength=length(indn);
ind=[];
ratio=(n-nn-1)/(nn-1);
for i=1:indLength;
    [k,l]=ind2sub([nn,nn],indn(i));
    k=k-1;
    l=l-1;
    if k==0 || l==0
        continue;
    end
    kb=floor((k-1)*ratio)+2;
    kt=floor(k*ratio)+1;
    lb=floor((l-1)*ratio)+2;
    lt=floor(l*ratio)+1;
    for kk=kb:kt
        for ll=lb:lt
            tmp=sub2ind([n-nn,n-nn], kk, ll);
            ind=[ind tmp];
        end
    end
end
% ind=[ind (n-nn)^2]; 

function  [p1]=PermutationTest(C,D,rep,option)
% This is an auxiliary function of the main function to calculate the p-values of
% all local tests of mcorr/dcorr/Mantel, the p-value of HHG in the
% permutation test.
if nargin<4
    option=1;  % Default option. Set to 1/2/3/4 for local mcorr/dcorr/Mantel, or HHG.
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
    dCor1=zeros(n,n,rep);  % Default option. Set to 1/2/3/4 for local mcorr/dcorr/Mantel, or HHG.
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
p1(1,:)=1;p1(:,1)=1; % Set the p-values of rank 0 to maximum.
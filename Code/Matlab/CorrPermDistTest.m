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
rl=50;
epsilon=0.5;

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
    rr=ceil(rep2/rep1);
    rr=max(rl,rr);
    ind1=splitPermutation(C,D,nn,rep1,rr,epsilon,option(1));
    ind2=splitPermutation(C,D,nn,rep1,rr,epsilon,option(2));
    ind3=splitPermutation(C,D,nn,rep1,rr,epsilon,option(3));
end
% % p1=median(p1All(ind1));
% % p2=median(p2All(ind2));
% % p3=median(p3All(ind3));
p1=mean(p1All(ind1));
p2=mean(p2All(ind2));
p3=mean(p3All(ind3));
% p1=min(p1All(ind1));
% p2=min(p2All(ind2));
% p3=min(p3All(ind3));

if p2<0.05 && p5>0.05
    length(ind2)
    p2
    p5
end
% if min(p1,p4)<0.05
%     length(ind1)
%     p1
%     p4
% end
% if min(p2,p5)<0.05
%     length(ind2)
%     p2
%     p5
% end
% Save the results
pre1='../../Data/';
filename=strcat(pre1,'CorrPermDistTestType',titlechar);
save(filename,'titlechar','rep2','option','p1All','p2All','p3All','p4','p5','p6','p7','p1','p2','p3','ind1','ind2','ind3');

function ind=splitPermutation(C,D,nn,rep1,rep2,epsilon,option)
ind=[];
if option==0;
    ind=1;
    return;
end
n=size(C,1);
ind=1:nn^2;
repeatInd=0;
dCorA=zeros(nn,nn,floor(n/nn)*rep1);
for r=1:rep1
    per=randperm(n);
    l1=length(ind);
    for rr=1:floor(n/nn)
        pa1=per((rr-1)*nn+1:rr*nn);
        [tmp]=PermutationTest(C(pa1,pa1),D(pa1,pa1),rep2,option);
        %         tmp=tmp-min(min(tmp));
        eps=min(min(min(tmp))+epsilon,0.8);
        ind=ind(tmp(ind)<eps);
        dCorA(:,:,(r-1)+rr)=tmp;
    end
    %     if length(ind)<round(nn);
    %         ind=verifyAdjacency(ind,nn);
    %     end
    if length(ind)==l1
        repeatInd=repeatInd+1;
    end
    if length(ind)~=l1
        repeatInd=0;
    end
    if repeatInd==round(rep1/6);
        break;
    end
    if length(ind)<nn;
        ind=verifyAdjacency(ind,nn);
        break;
    end
    %     if length(ind)<nn/10;%round(0.1*nn);
    %         ind=[];
    %         break;
    %     end
end
if isempty(ind);
    ind=n^2;
else
    ind=ratioScale(ind,nn,n);
end

function ind=verifyAdjacency(ind,nn)

% if length(ind)<nn
ll=length(ind);
row=zeros(ll,1);
col=zeros(ll,1);
for i=1:ll
    [k,l]=ind2sub([nn,nn],ind(i));
    row(i)=k;
    col(i)=l;
end
[km,ck]=mode(row);
[lm,cl]=mode(col);
tt=false;
if ck>ll/1.5 && ck>cl
    pos=(row==km);
    ind=ind(pos);
    tt=true;
end
if cl>ll/1.5
    pos=(col==lm);
    ind=ind(pos);    
    tt=true;
end
% length(ind)
if length(ind)<nn/5 || tt==false
    ind=[];
end

function ind=splitPermutation2(C,D,nn,rep1,rep2,epsilon,option)
ind=[];
if option==0;
    ind=1;
    return;
end
n=size(C,1);
ind=1:nn^2;
dCorA=zeros(floor(n/nn)*rep1,nn,nn);
for r=1:rep1
    per=randperm(n);
    for rr=1:floor(n/nn)
        pa1=per((rr-1)*nn+1:rr*nn);
        [tmp]=PermutationTest(C(pa1,pa1),D(pa1,pa1),rep2,option);
%         tmp=tmp-min(min(tmp));
%         ind=ind(tmp(ind)<epsilon);
        dCorA(floor(n/nn)*(r-1)+rr,:,:)=tmp;
    end
end

bench=dCorA(:,end,end);
h=zeros(nn,nn);
unif=0:0.01:1;
if kstest2(unif,bench,'Tail','larger')==1;
    bench=unif;
end
for i=1:nn
    for j=1:nn
       [~,h(i,j)]=kstest2(dCorA(:,i,j),bench,'Tail','larger');
    end
end
ind=find(h<0.05);

pos=[];
tmp=1;
for i=1:length(ind)
    [k,l]=ind2sub([nn,nn],ind(i));
    tmp=min(tmp,max(dCorA(:,k,l)));
    if max(dCorA(:,k,l))<epsilon
        pos=[pos i];
    end
end
tmp
if isempty(pos)
    ind=n^2;
else
    ind=ind(pos);
    ind=ratioScale(ind,nn,n);
end
% end

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
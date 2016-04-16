function [p1,p2,p3, p4,p5,p6,p7,p1All,p2All,p3All]=CorrPermDistTest(C,D,rep,titlechar,option)
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
if nargin<4
    titlechar='RealData';
end
if nargin<5
    option=[1,2,3,4]; % Default option. Setting any to 0 to disable the calculation of MGC by dcorr/mcorr/Mantel, global dcorr/mcorr/Mantel, HHG, in order; set the first three
end

% Global Correlations
[p1All]=PermutationTest(C,D,rep,option(1));
[p2All]=PermutationTest(C,D,rep,option(2));
[p3All]=PermutationTest(C,D,rep,option(3));
[p1,ind1]=verify(p1All);[p2,ind2]=verify(p2All);[p3,ind3]=verify(p3All);
p4=p1All(end);p5=p2All(end);p6=p3All(end);
[p7]=PermutationTest(C,D,rep,option(4));

if p2<0.05 || p5<0.05
    p2
    p5
end
pre1='../../Data/';
filename=strcat(pre1,'CorrPermDistTestType',titlechar);
save(filename,'titlechar','rep','option','p1All','p2All','p3All','p4','p5','p6','p7','p1','p2','p3','ind1','ind2','ind3');

function [p1,ind]=verify(p1All)
alpha=0.05;
thres=0.80;
if p1All==1;
    p1=1;
    ind=[];
    return;
end
n=size(p1All,1);
power1=(p1All<alpha);
pCol=mean(power1(2:end,1:end),1);
pRow=mean(power1(1:end,2:end),2);
% figure
% plot(1:n,pCol,'bo:',1:n,pRow,'r+:')

if length(find(pCol>thres))+length(find(pRow>thres))<n/20 %&& length(find(pRow>thres2))<n/20
%     figure
%     plot(1:n,pCol,'bo:',1:n,pRow,'r+:')

    ind=[];
    p1=1;
%     ind=n^2;
%     p1=p1All(ind);
    return;
end

% figure
% imagesc(pMean)
% colorbar();
% figure
% plot(1:nn,pCol,'bo:',1:nn,pRow,'r+:')
if max(pRow)>max(pCol)
    k=find(pRow==max(pRow),1,'last');
    l=find(p1All(k,:)==min(p1All(k,:)),1,'last');
else
    l=find(pCol==max(pCol),1,'last');
    k=find(p1All(:,l)==min(p1All(:,l)),1,'last');
end
ind=sub2ind([n,n],k,l);
p1=p1All(ind);

function [p1,ind]=verify3(p1All,thres)
if nargin<2
    thres=0.05;
end
if p1All==1;
    p1=1;
    ind=[];
    return;
end
p1=0;
n=size(p1All,1);
pCol=mean(p1All(2:end,1:end),1);
pRow=mean(p1All(1:end,2:end),2);

if min(pCol)>thres && min(pRow)>thres
    ind=[];
    p1=1;
    return;
end

if min(pRow)>min(pCol)
    p1All=p1All';
    pCol=mean(p1All(2:end,1:end),1);
    pRow=mean(p1All(1:end,2:end),2);
end

pos=find(pRow<thres);
ind=find(pRow(pos)==min(pRow(pos)),1,'last');
if isempty(ind)
    return;
end
k=pos(ind);
l=find(p1All(k,:)==min(p1All(k,:)),1,'last');
% figure
% plot(1:nn,pCol,'bo:',1:nn,pRow,'r+:')
% figure
% imagesc(pMean);
% colorbar();
% length(pos)
if length(pos)<n/20;
    ind=[];
    p1=1;
    return;
else
    ind=sub2ind([n,n],k,l);
    p1=p1All(ind);
end

function pos=verifyAdjacentScales(pos)
pos2=find(diff(pos)>=0);
pos1=verifyAuxi(pos2);

pos2=find(diff(pos)<=0);
pos2=verifyAuxi(pos2);

if length(pos1)>length(pos2)
    pos2=pos1;
end
if isempty(pos2)
    pos=pos(end);
else
    pos=pos(pos2);
end

function pos2=verifyAuxi(pos2)
pos2=find(diff(pos2)==1);
if isempty(pos2)
    pos2=[];
    return;
end
posStart=1;
posStartTmp=1;
maxT=1;
maxTmp=1;
for i=2:length(pos2)
    if pos2(i)==pos2(i-1)+1
        maxTmp=maxTmp+1;
        if maxTmp>=maxT
            maxT=maxTmp;
            posStart=posStartTmp;
        end
    else
        posStartTmp=i;
        maxTmp=1;
    end
end
pos2=pos2(posStart):pos2(posStart+maxTmp-1)+1;

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
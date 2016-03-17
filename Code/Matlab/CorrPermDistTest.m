function [p1, p2, p3, p4,p5,p6,p7,neighborhoods]=CorrPermDistTest(C,D,rep1,rep2,titlechar,alpha,option,neighborhoods)
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
    alpha=0.05; % Default type 1 error level
end
if nargin<7
    option=[1,1,1,1,1,1,1]; % Default option. Setting any to 0 to disable the calculation of MGC by mcorr/dcorr/Mantel, global mcorr/dcorr/Mantel, HHG, in order; set the first three
end
if nargin<8
    neighborhoods=ones(3,1);
end

% % Run the independence test to first estimate the optimal scale of MGC
% if rep1~=0 && mean(neighborhoods-1)<=0
%     [p1,p2,p3,n1,n2,n3]=IndependenceTest(C,D,rep1,option,alpha);
%     neighborhoods=[n1;n2;n3];
%     p4=0;p5=0;p6=0;p7=0;
% end

% Run the permutation test to return the p-values
% if rep2~=0
    [p1All, p2All, p3All, p7,n1,n2,n3]=PermutationTest(C,D,rep1,rep2,option,alpha);
    % From the p-values of all local tests,  get the p-values of MGC based on the optimal neighborhood estimation, and the p-values of the respective global test
    if mean(neighborhoods-1)<=0
        neighborhoods=[n1;n2;n3];
    end
%     if mean(neighborhoods-1)>0
        p1=p1All(neighborhoods(1));
        p2=p2All(neighborhoods(2));
        p3=p3All(neighborhoods(3));
%     else
%         p1=min(min(p1All));p2=min(min(p2All));p3=min(min(p3All));
%     end
    p4=p1All(end);p5=p2All(end);p6=p3All(end);
    % Save the results
    pre1='../../Data/';
    filename=strcat(pre1,'CorrPermDistTestType',titlechar);
    save(filename,'titlechar','p1','p2','p3','p4','p5','p6','p7','neighborhoods','rep1','rep2','alpha','option','p1All','p2All','p3All');
% end

function  [power1, power2, power3,n1,n2,n3]=IndependenceTest2(C,D,rep,option,alpha)
% This is an auxiliary function of the main function to calculate the powers of
% all local tests of mcorr/dcorr/Mantel, in order to locate the optimal
% scale of MGC.
%
% By kernel density estimation on the upper-diagonal entries of the distance
% matrices, it independently generates paired distance matrices that
% distribute similarly as the given pair (C,D). Then the test statistics under the null and the
% alternative can be calculated to estimate the testing power.
n=size(C,1);
% Test statiscs under the null and the alternative
dCor1N=zeros(n,n,rep);dCor2N=zeros(n,n,rep);dCor3N=zeros(n,n,rep);
dCor1A=zeros(n,n,rep);dCor2A=zeros(n,n,rep);dCor3A=zeros(n,n,rep);
% Powers for all local tests of mcorr/dcorr/Mantel
power1=zeros(n,n);power2=zeros(n,n);power3=zeros(n,n);

% % Kernal density estimation on the lower-diagonal entries of the distance matrices
CU=C;DU=D;
H=eye(n)-(1/n)*ones(n,n);
CU=H*C*H;DU=H*D*H;
meanC=C-CU; meanD=D-DU;
ind=zeros(n*(n-1)/2,1);
for i=1:n;
    for j=1:i-1;
        ind((i-1)*(i-2)/2+j)=n*(j-1)+i;
    end
end
CU=CU(ind);DU=DU(ind);
mC=mean(CU);mD=mean(DU);
varC=std(CU)^2;varD=std(DU)^2;
[bandwidth]=kde2d([CU, DU],256);
% If the kernal density estimation fails, use random resampling instead
if isnan(bandwidth)
    %disp('Kernal density estimation failed; use random sampling instead');
    bandwidth=[0;0];
end
sz1=size(CU,1);
sz2=n*(n-1)/2;

for r=1:rep
    % Generate an independent pair of the upper diagonal distance entries
%     aa=rand(sz2,1);bb=randn(sz2,1);
%     CAU = mC+(CU(ceil(sz1*aa)) -mC + bandwidth(1)*bb)/sqrt(1+bandwidth(1)^2/varC);
%     DAU = mD+(DU(ceil(sz1*aa)) -mD + bandwidth(2)*bb)/sqrt(1+bandwidth(2)^2/varD);
%     % Re-form the symmetric distance matrices. The triangle inequality may
%     % not hold though.
%     CA=zeros(n,n);DA=zeros(n,n);
%     for i=2:n;
%         for j=1:i-1;
%             indd=(i-1)*(i-2)/2;
%             CA(i,j)=CAU(indd+j);
%             DA(i,j)=DAU(indd+j);
%             CA(j,i)=CA(i,j);
%             DA(j,i)=DA(i,j);
%         end
%     end
%     CA=CA+meanC;DA=DA+meanD;
    noise1=squareform(pdist(rand(n,1)));
    noise1=H*noise1*H;
    CA=C+1*noise1;
    DA=D;%+10*noise1;
    % Calculate the test statistics under the alternative
    disRank1=[disToRanks(CA) disToRanks(DA)];
    if option(1)~=0
        dCor1A(:,:,r)=LocalGraphCorr(CA,DA,1,disRank1);
    end
    if option(2)~=0
        dCor2A(:,:,r)=LocalGraphCorr(CA,DA,2,disRank1);
    end
    if option(3)~=0
        dCor3A(:,:,r)=LocalGraphCorr(CA,DA,3,disRank1);
    end
    
    % Generate another independent second distance matrix
%     per1=randsample(n,n,true);
%     per2=randsample(n,n,true);
    per1=randperm(n);
    per2=randperm(n);
    DN=DA(per1,per1);
    CN=CA(per2,per2);
    % Calculate the test statistics under the null
    disRank1=[disToRanks(CN) disToRanks(DN)];
    if option(1)~=0
        dCor1N(:,:,r)=LocalGraphCorr(CN,DN,1,disRank1);
    end
    if option(2)~=0
        dCor2N(:,:,r)=LocalGraphCorr(CN,DN,2,disRank1);
    end
    if option(3)~=0
        dCor3N(:,:,r)=LocalGraphCorr(CN,DN,3,disRank1);
    end
end

% % For each local test, estimate the critical value from the test statistics under the null,
% % then estimate the power from the test statistics under the alternative.
for i=1:n
    for j=1:n;
        dCorT=sort(dCor1N(i,j,:),'descend');
        cut1=dCorT(ceil(rep*alpha));
        power1(i,j)=mean(dCor1A(i,j,:)>cut1);
        
        dCorT=sort(dCor2N(i,j,:),'descend');
        cut2=dCorT(ceil(rep*alpha));
        power2(i,j)=mean(dCor2A(i,j,:)>cut2);
        
        dCorT=sort(dCor3N(i,j,:),'descend');
        cut3=dCorT(ceil(rep*alpha));
        power3(i,j)=mean(dCor3A(i,j,:)>cut3);
    end
end
% Set the powers of all local tests at rank 0 to 0
power1(1,:)=0;power1(:,1)=0;power2(1,:)=0;power2(:,1)=0;power3(1,:)=0;power3(:,1)=0;
n1=maxNeighbors(power1,dCor1N,dCor1A);
n2=maxNeighbors(power2,dCor2N,dCor2A);
n3=maxNeighbors(power3,dCor3N,dCor3A);
save('tmpPerm.mat','dCor1N','dCor2N','dCor3N','dCor1A','dCor2A','dCor3A','power1','power2','power3','n1','n2','n3');

function  [p1, p2, p3,p4,n1,n2,n3]=PermutationTest(C,D,rep1,rep2,option,alpha)
if nargin<5
    option=[1,1,1,1];  % Default option. Setting each entry to 0 to disable the calculation of local mcorr/dcorr/Mantel, or HHG.
end
n=size(C,1);
epsilon=0.001;

% P-values of all local tests of mcorr/dcorr/Mantel
p1=zeros(n,n); p2=zeros(n,n);p3=zeros(n,n);p4=0;
dCor1=zeros(n,n,rep2);dCor2=zeros(n,n,rep2);dCor3=zeros(n,n,rep2);

% Calculate the observed test statistics for the given data sets
disRankC=disToRanks(C);
disRankP=disToRanks(D);
disRank=[disRankC disRankP];
if option(1)~=0
    cut1=LocalGraphCorr(C,D,1,disRank);
end
if option(2)~=0
    cut2=LocalGraphCorr(C,D,2,disRank);
end
if option(3)~=0
    cut3=LocalGraphCorr(C,D,3,disRank);
end
if option(4)~=0
    cut4=HHG(C,D);
end

% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep2
    % Use random permutations;
    per=randperm(n);
    DN=D(per,per);
    CN=C;
    disRank=[disToRanks(CN) disToRanks(DN)];
    if option(1)~=0
        tmp=LocalGraphCorr(CN,DN,1,disRank);
        p1=p1+(tmp>=cut1)/rep2;
        dCor1(:,:,r)=tmp;
    end
    if option(2)~=0
        tmp=LocalGraphCorr(CN,DN,2,disRank);
        p2=p2+(tmp>=cut2)/rep2;
        dCor2(:,:,r)=tmp;
    end
    if option(3)~=0
        tmp=LocalGraphCorr(CN,DN,3,disRank);
        p3=p3+(tmp>=cut3)/rep2;
        dCor3(:,:,r)=tmp;
    end
    if option(4)~=0
        dCor4=HHG(CN,DN);
        p4=p4+(dCor4>=cut4)/rep2;
    end
end

n1=maxNeighbors(1-p1);
n2=maxNeighbors(1-p2);
n3=maxNeighbors(1-p3);

pp1=zeros(rep1,1);
pp2=zeros(rep1,1);
pp3=zeros(rep1,1);
for r=1:rep1
    % Use random permutations;
    if option(1)~=0
        tmp=0;
        tmpC=dCor1(:,:,r);
        for rr=1:rep1
             tmp=tmp+(dCor1(:,:,rr)>=tmpC)/rep1;
        end
        pp1(r)=min(min(tmp(2:n,2:n)))-1/rep1;
        % pp1(r)=tmp(n1);
    end
    if option(2)~=0
        tmp=0;
        tmpC=dCor2(:,:,r);
        for rr=1:rep1
             tmp=tmp+(dCor2(:,:,rr)>=tmpC)/rep1;
        end
        pp2(r)=min(min(tmp(2:n,2:n)))-1/rep1;
    end
    if option(3)~=0
        tmp=0;
        tmpC=dCor3(:,:,r);
        for rr=1:rep1
             tmp=tmp+(dCor3(:,:,rr)>=tmpC)/rep1;
        end
        pp3(r)=min(min(tmp(2:n,2:n)))-1/rep1;
    end
end

pp1=sort(pp1,'ascend');
pp2=sort(pp2,'ascend');
pp3=sort(pp3,'ascend');
%p1(n1)
%pp1(floor(alpha*rep1))
p1(n1)=p1(n1)/(pp1(floor(alpha*rep1))+epsilon)*alpha;
p2(n2)=p2(n2)/(pp2(floor(alpha*rep1))+epsilon)*alpha;
p3(n3)=p3(n3)/(pp3(floor(alpha*rep1))+epsilon)*alpha;
% if p1(n1)<pp1(ceil(alpha*rep2));
%     n1=n^2;
% end
% if p2(n1)<pp2(ceil(alpha*rep2));
%     n2=n^2;
% end
% if p3(n1)<pp3(ceil(alpha*rep2));
%     n3=n^2;
% end
% Output the p-value

%save('tmpPerm.mat','dCor1','dCor2','dCor3','p1','p2','p3','n1','n2','n3');

% function  [p1, p2, p3, p4]=PermutationTest(C,D,rep,option)
% % This is an auxiliary function of the main function to calculate the p-values of
% % all local tests of mcorr/dcorr/Mantel, the p-value of HHG in the
% % permutation test.
% if nargin<5
%     option=[1,1,1,1];  % Default option. Setting each entry to 0 to disable the calculation of local mcorr/dcorr/Mantel, or HHG.
% end
% n=size(C,1);
% 
% % P-values of all local tests of mcorr/dcorr/Mantel
% p1=zeros(n,n); p2=zeros(n,n);p3=zeros(n,n);
% p4=0; % P-values for HHG
% 
% % Calculate the observed test statistics for the given data sets
% disRankC=disToRanks(C);
% disRankP=disToRanks(D);
% disRank=[disRankC disRankP];
% if option(1)~=0
%     cut1=LocalGraphCorr(C,D,1,disRank);
% end
% if option(2)~=0
%     cut2=LocalGraphCorr(C,D,2,disRank);
% end
% if option(3)~=0
%     cut3=LocalGraphCorr(C,D,3,disRank);
% end
% if option(4)~=0
%     cut4=HHG(C,D);
% end
% 
% % Now Permute the second dataset for rep times, and calculate the p-values
% for r=1:rep
%     % Use random permutations;
%     per=randperm(n);
%     DN=D(per,per);
%     CN=C;
%     disRank=[disToRanks(CN) disToRanks(DN)];
%     if option(1)~=0
%         dCor1=LocalGraphCorr(CN,DN,1,disRank);
%         p1=p1+(dCor1<cut1)/rep;
%     end
%     if option(2)~=0
%         dCor2=LocalGraphCorr(CN,DN,2,disRank);
%         p2=p2+(dCor2<cut2)/rep;
%     end
%     if option(3)~=0
%         dCor3=LocalGraphCorr(CN,DN,3,disRank);
%         p3=p3+(dCor3<cut3)/rep;
%     end
%     if option(4)~=0
%         dCor4=HHG(CN,DN);
%         p4=p4+(dCor4<cut4)/rep;
%     end
% end
% 
% % Output the p-value
% p1=1-p1;p2=1-p2;p3=1-p3;p4=1-p4;
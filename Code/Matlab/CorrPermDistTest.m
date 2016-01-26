function [p1, p2, p3, p4,p5,p6,p7,neighborhoods]=CorrPermDistTest(C,P,rep1,rep2,titlechar,neighborhoods,allP,alpha,option)
% Author: Cencheng Shen
% Permutation Tests for identifying dependency.
% The output are the p-values of MGC by mcorr/dcorr/Mantel, and global mcorr/dcorr/Mantel/HHG.

% Parameters:
% type should be a n*2n matrix, a concatenation of two distance matrices,
% rep1 specifies the number of bootstrap samples for optimal scale estimation
% in MGC; if 0, the estimation step is skipped.
% rep2 specifies the number of random permutations to use for the permutation test,
% set allP to non-zero will use all permutations instead,
% alpha specifies the type 1 error level,
% option specifies whether each test statistic is calculated or not.
if nargin<5
    titlechar='RealData';
end
if nargin<6
    neighborhoods=zeros(3,1); % Estimated optimal neighborhoods at each sample size. At 0, MGC of all scales are calculated
end
if nargin<7
    allP=0; % If set to other value, will override rep and use all permutations; unfeasible for n large.
end
if nargin<8
    alpha=0.05; % Default type 1 error level
end
if nargin<9
    option=[1,1,1,1,1,1,1]; % Default option. Setting any to 0 to disable the calculation of MGC by mcorr/dcorr/Mantel, global mcorr/dcorr/Mantel, HHG, in order; set the first three
end

% Run the independence test to first estimate the optimal scale of MGC
if rep1~=0
    [p1,p2,p3]=IndependenceTest(C,P,rep1,alpha);
    if neighborhoods(1)==0
        neighborhoods(1)=verifyNeighbors(1-p1);
    end
    if neighborhoods(2)==0
        neighborhoods(2)=verifyNeighbors(1-p2);
    end
    if neighborhoods(3)==0
        neighborhoods(3)=verifyNeighbors(1-p3);
    end
    p4=0;p5=0;p6=0;p7=0;
end

% Return the permutation test to return the p-values
if rep2~=0
    [p1All, p2All, p3All, p7]=PermutationTest(C,P,rep2,allP,option);
    if rep1==0
        if neighborhoods(1)==0
            neighborhoods(1)=verifyNeighbors(p1All);
        end
        if neighborhoods(2)==0
            neighborhoods(2)=verifyNeighbors(p2All);
        end
        if neighborhoods(3)==0
            neighborhoods(3)=verifyNeighbors(p3All);
        end
    end
    % From the p-values of all local tests,  get the p-values of MGC based on the optimal neighborhood estimation, and the p-values of the respective global test
    p1=p1All(neighborhoods(1));p4=p1All(end);
    p2=p2All(neighborhoods(2));p5=p2All(end);
    p3=p3All(neighborhoods(3));p6=p3All(end);
    % Save the results
    pre1='../../Data/';
    filename=strcat(pre1,'CorrPermDistTestType',titlechar);
    save(filename,'titlechar','p1','p2','p3','p4','p5','p6','p7','neighborhoods','rep1','rep2','allP','alpha','option','p1All','p2All','p3All');
end

% Smooth Bootstrap
function  [power1, power2, power3]=IndependenceTest(C,P,rep,alpha)
% Author: Cencheng Shen
% This is an auxiliary function of the main function to calculate the powers of
% all local tests of mcorr/dcorr/Mantel, in order to locate the optimal
% scale of MGC.
%
% It generates dependent and independent data from smooth bootstrap samples,
% calculate the test statistics under the null and the alternative,
% then estimate the testing power of each method.
n=size(C,1);
% Test statiscs under the null and the alternative
dCor1N=zeros(n,n,rep);dCor2N=zeros(n,n,rep);dCor3N=zeros(n,n,rep);
dCor1A=zeros(n,n,rep);dCor2A=zeros(n,n,rep);dCor3A=zeros(n,n,rep);
% Powers for all local tests of mcorr/dcorr/Mantel
power1=zeros(n,n);power2=zeros(n,n);power3=zeros(n,n);
ratio=1/n;
    
for r=1:rep
    % Random sampling with replacement
    per=randsample(n,n,true);
    noise1=random('norm',0,1,n,1);
    noise1=squareform(pdist(noise1));
    noise1=noise1/norm(noise1,'fro')/norm(C,'fro');
    Ca=C(per,per);
    %Ca=sqrt(Ca.^2+noise1.^2*ratio+Ca.*noise1*2*sqrt(ratio));
    Ca=Ca+noise1*sqrt(ratio);
%     noise2=noise1/norm(noise1,'fro')/norm(P,'fro');
%     P=sqrt(P.^2+noise2.^2*ratio+P.*noise2*2*sqrt(ratio));
    %Ca=C(per,per);
%     noise1=random('norm',0,1,n,1);
%     noise1=squareform(pdist(noise1));
%     noise1=noise1/norm(noise1,'fro')/norm(P,'fro');
    Pa=P(per,per); 
%     Pa=sqrt(Pa.^2+noise1.^2*ratio+Pa.*noise1*2*sqrt(ratio));
    disRank=[disToRanks(Ca) disToRanks(Pa)];
    % Calculate the test statistics under the alternative
    dCor1A(:,:,r)=LocalGraphCorr(Ca,Pa,1,disRank);
    dCor2A(:,:,r)=LocalGraphCorr(Ca,Pa,2,disRank);
    dCor3A(:,:,r)=LocalGraphCorr(Ca,Pa,3,disRank);
    
    % A different random sampling
    perN=randsample(n,n,true);
%     noise1=random('norm',0,1,n,1);
%     noise1=squareform(pdist(noise1));
%     noise2=noise1/norm(noise1,'fro')/norm(P,'fro');
%     P=sqrt(P.^2+noise2.^2/n+P.*noise2*2/sqrt(n));
    Pa=Pa(perN,perN);
    disRank=[disToRanks(Ca) disToRanks(Pa)];
    % Calculate the test statistics under the null
    dCor1N(:,:,r)=LocalGraphCorr(Ca,Pa,1,disRank);
    dCor2N(:,:,r)=LocalGraphCorr(Ca,Pa,2,disRank);
    dCor3N(:,:,r)=LocalGraphCorr(Ca,Pa,3,disRank);
end

% For each local test, estimate the critical value from the test statistics under the null,
% then estimate the power from the test statistics under the alternative.
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

% Usual Bootstrap
function  [power1, power2, power3]=IndependenceTest2(C,P,rep,alpha)
% Author: Cencheng Shen
% This is an auxiliary function of the main function to calculate the powers of
% all local tests of mcorr/dcorr/Mantel, in order to locate the optimal
% scale of MGC.
%
% It generates dependent and independent data from bootstrap samples,
% calculate the test statistics under the null and the alternative,
% then estimate the testing power of each method.
n=size(C,1);
% Test statiscs under the null and the alternative
dCor1N=zeros(n,n,rep);dCor2N=zeros(n,n,rep);dCor3N=zeros(n,n,rep);
dCor1A=zeros(n,n,rep);dCor2A=zeros(n,n,rep);dCor3A=zeros(n,n,rep);
% Powers for all local tests of mcorr/dcorr/Mantel
power1=zeros(n,n);power2=zeros(n,n);power3=zeros(n,n);
for r=1:rep
    % Random sampling with replacement
    per=randsample(n,n,true);
    Ca=C(per,per);
    Pa=P(per,per);
    disRank=[disToRanks(Ca) disToRanks(Pa)];
    % Calculate the test statistics under the alternative
    dCor1A(:,:,r)=LocalGraphCorr(Ca,Pa,1,disRank);
    dCor2A(:,:,r)=LocalGraphCorr(Ca,Pa,2,disRank);
    dCor3A(:,:,r)=LocalGraphCorr(Ca,Pa,3,disRank);
    
    % A different random sampling
    perN=randsample(n,n,true);
    Pa=P(perN,perN);
    disRank=[disToRanks(Ca) disToRanks(Pa)];
    % Calculate the test statistics under the null
    dCor1N(:,:,r)=LocalGraphCorr(Ca,Pa,1,disRank);
    dCor2N(:,:,r)=LocalGraphCorr(Ca,Pa,2,disRank);
    dCor3N(:,:,r)=LocalGraphCorr(Ca,Pa,3,disRank);
end

% For each local test, estimate the critical value from the test statistics under the null,
% then estimate the power from the test statistics under the alternative.
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
%
% % Set the powers of all local tests at rank 0 to alpha
% power1(1,:)=0;power1(:,1)=0;power2(1,:)=0;power2(:,1)=0;power3(1,:)=0;power3(:,1)=0;

function  [p1, p2, p3, p4]=PermutationTest(C,P,rep,allP,option)
% Author: Cencheng Shen
% This is an auxiliary function of the main function to calculate the p-values of
% all local tests of mcorr/dcorr/Mantel, the p-value of HHG in the
% permutation test.
if nargin<5
    option=[1,1,1,1];  % Default option. Setting each entry to 0 to disable the calculation of local mcorr/dcorr/Mantel, or HHG.
end
n=size(C,1);
if allP~=0
    PAll=perms(1:n);
    rep=size(PAll,1);
end

% P-values of all local tests of mcorr/dcorr/Mantel
p1=zeros(n,n); p2=zeros(n,n);p3=zeros(n,n);
p4=0; % P-values for HHG

% Calculate the observed test statistics for the given data sets
disRankC=disToRanks(C);
disRankP=disToRanks(P);
disRank=[disRankC disRankP];
if option(1)~=0
    cut1=LocalGraphCorr(C,P,1,disRank);
end
if option(2)~=0
    cut2=LocalGraphCorr(C,P,2,disRank);
end
if option(3)~=0
    cut3=LocalGraphCorr(C,P,3,disRank);
end
if option(4)~=0
    cut4=HHG(C,P);
end

% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % Use random permutations; if allP is not 0, use all possible permutations
    per=randperm(n);
    if allP~=0
        per=PAll(r,:);
    end
    Pa=P(per,per);
    disRank=[disRankC disRankP(per, per)];
    if option(1)~=0
        dCor1=LocalGraphCorr(C,Pa,1,disRank);
        p1=p1+(dCor1<cut1)/rep;
    end
    if option(2)~=0
        dCor2=LocalGraphCorr(C,Pa,2,disRank);
        p2=p2+(dCor2<cut2)/rep;
    end
    if option(3)~=0
        dCor3=LocalGraphCorr(C,Pa,3,disRank);
        p3=p3+(dCor3<cut3)/rep;
    end
    if option(4)~=0
        dCor4=HHG(C,Pa);
        p4=p4+(dCor4<cut4)/rep;
    end
end

% Output the p-value
p1=1-p1;p2=1-p2;p3=1-p3;p4=1-p4;

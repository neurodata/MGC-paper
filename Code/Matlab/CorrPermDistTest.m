function [p1, p2, p3, p4,p5,p6,p7]=CorrPermDistTest(type,rep1,rep2,titlechar, allP, option)
% Author: Cencheng Shen
% Permutation Tests for identifying dependency.
% The output are the p-values of MGC by mcorr/dcorr/Mantel, and global mcorr/dcorr/Mantel/HHG.

% Parameters:
% type should be a n*2n matrix, a concatenation of two distance matrices,
% rep1 specifies the number of bootstrap samples for optimal scale estimation
% in MGC; if 0, the estimation step is skipped and the p-values of all local
% tests are returned instead.
% rep2 specifies the number of random permutations to use for the permutation test,
% set allP to non-zero will use all permutations instead,
% option specifies whether each test statistic is calculated or not.
if nargin<4
    titlechar='Real Data';
end
if nargin<5
    allP=0; % If set to other value, will override rep and use all permutations; unfeasible for n large.
end
if nargin<6
    option=[1,1,1,1,1,1,1]; % Default option. Setting any to 0 to disable the calculation of MGC by mcorr/dcorr/Mantel, global mcorr/dcorr/Mantel, HHG, in order.
end
n=size(type,1);
C=type(:, 1:n);
P=type(:, n+1:2*n);

neighborhoods=zeros(3,1); % Estimated optimal neighborhoods at each sample size. At 0, MGC of all scales are calculated
    
% Run the independence test to first estimate the optimal scale of MGC
if rep1~=0
    [power1,power2,power3]=IndependenceTest(C,P,rep1);
    neighbor1=verifyNeighbors(1-power1);
    neighbor2=verifyNeighbors(1-power2);
    neighbor3=verifyNeighbors(1-power3);
    neighborhoods=[neighbor1(end),neighbor2(end),neighbor3(end)];
end
% Return the permutation test to return the p-values
[p1, p2, p3, p4,p5,p6,p7]=PermutationTest(C,P,rep2,allP,option,neighborhoods);
if rep1==0
    neighbor1=verifyNeighbors(p1);
    neighbor2=verifyNeighbors(p2);
    neighbor3=verifyNeighbors(p3);
    neighborhoods=[neighbor1(end),neighbor2(end),neighbor3(end)];
end
% Save the results
if rep1==0
    tmpC='All'; 
else
    tmpC='';
end
filename=strcat('CorrPermDistTestType',titlechar,tmpC);
save(filename,'titlechar','p1','p2','p3','p4','p5','p6','p7','neighborhoods','type','n','rep1','rep2','allP','option');

function  [p1, p2, p3, p4,p5,p6,p7]=PermutationTest(C,P,rep,allP,option,neighborhoods)
% Author: Cencheng Shen
% This is an auxiliary function of the main test.
%
% It returs the p-value of MGC by mcorr/dcorr/Mantel, and global
% mcorr/dcorr/Mantel/HHG
if nargin<5
    option=[1,1,1,1,1,1,1];  % Default option. Setting any to 0 to disable the calculation of MGC by mcorr/dcorr/Mantel, global mcorr/dcorr/Mantel, HHG, in order.
end
if nargin<6 || length(neighborhoods)<3
    neighborhoods=zeros(3,1); % Default neighborhood to computes MGC of all scales
end
if allP~=0 
    PAll=perms(1:n);
    rep=size(PAll,1);
end

n=size(C,1);
if norm(neighborhoods,'fro')==0
    indD=[n,n];
else
    indD=1;
end
p1=zeros(indD); p2=zeros(indD);p3=zeros(indD); % P-value of MGC by mcorr/dcorr/Mantel
p4=0;p5=0;p6=0;p7=0; % P-values for global mcorr/dcorr/Mantel/HHG

% Calculate the observed test statistics for the given data sets
disRankC=disToRanks(C);
disRankP=disToRanks(P);
disRank=[disRankC disRankP];
if option(1)~=0
    cut1=LocalGraphCorr(C,P,1,neighborhoods(1),disRank);
end
if option(2)~=0
    cut2=LocalGraphCorr(C,P,2,neighborhoods(2),disRank);
end
if option(3)~=0
    cut3=LocalGraphCorr(C,P,3,neighborhoods(3),disRank);
end
if option(4)~=0
    cut4=LocalGraphCorr(C,P,1,n^2,disRank);
end
if option(5)~=0
    cut5=LocalGraphCorr(C,P,2,n^2,disRank);
end
if option(6)~=0
    cut6=LocalGraphCorr(C,P,3,n^2,disRank);
end
if option(7)~=0
    cut7=HHG(C,P);
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
        dCor1=LocalGraphCorr(C,Pa,1,neighborhoods(1),disRank);
        p1=p1+(dCor1<cut1)/rep;
    end
    if option(2)~=0
        dCor2=LocalGraphCorr(C,Pa,2,neighborhoods(2),disRank);
        p2=p2+(dCor2<cut2)/rep;
    end
    if option(3)~=0
        dCor3=LocalGraphCorr(C,Pa,3,neighborhoods(3),disRank);
        p3=p3+(dCor3<cut3)/rep;
    end
    if option(4)~=0
        dCor4=LocalGraphCorr(C,Pa,1,n^2,disRank);
        p4=p4+(dCor4<cut4)/rep;
    end
    if option(5)~=0
        dCor5=LocalGraphCorr(C,Pa,2,n^2,disRank);
        p5=p5+(dCor5<cut5)/rep;
    end
    if option(6)~=0
        dCor6=LocalGraphCorr(C,Pa,3,n^2,disRank);
        p6=p6+(dCor6<cut6)/rep;
    end
    if option(7)~=0
        dCor7=HHG(C,Pa);
        p7=p7+(dCor7<cut7)/rep;
    end
end

% Output the p-value
p1=1-p1;
p2=1-p2;
p3=1-p3;
p4=1-p4;
p5=1-p5;
p6=1-p6;
p7=1-p7;

function  [power1, power2, power3]=IndependenceTest(C,P,rep)
% Author: Cencheng Shen
% This is an auxiliary function of the main test.
%
% It generates dependent and independent data from resampling, 
% calculate the test statistics under the null and the alternative, 
% then estimate the testing power of each method.
% 
% Note that this is only for MGC to find the optimal scale, thus it is
% only for local tests and no global test

n=size(C,1);
alpha=0.05; % type 1 error level
% Test statiscs under the null and the alternative
dCor1N=zeros(n,n,rep);dCor2N=zeros(n,n,rep);dCor3N=zeros(n,n,rep);
dCor1A=zeros(n,n,rep);dCor2A=zeros(n,n,rep);dCor3A=zeros(n,n,rep);
% Powers for MGC by mcorr/dcorr/Mantel
power1=zeros(n,n);power2=zeros(n,n);power3=zeros(n,n);

for r=1:rep
    % Random sampling with replacement
    per=randsample(n,n,true);
    Ca=C(per,per);
    Pa=P(per,per);
    disRank=[disToRanks(Ca) disToRanks(Pa)];
    % Calculate the test statistics under the alternative
    dCor1A(:,:,r)=LocalGraphCorr(Ca,Pa,1,0,disRank);
    dCor2A(:,:,r)=LocalGraphCorr(Ca,Pa,2,0,disRank);
    dCor3A(:,:,r)=LocalGraphCorr(Ca,Pa,3,0,disRank);
    
    % A different random sampling
    perN=randsample(n,n,true);
    Pa=P(perN,perN);
    disRank=[disToRanks(Ca) disToRanks(Pa)];
    % Calculate the test statistics under the null
    dCor1N(:,:,r)=LocalGraphCorr(Ca,Pa,1,0,disRank);
    dCor2N(:,:,r)=LocalGraphCorr(Ca,Pa,2,0,disRank);
    dCor3N(:,:,r)=LocalGraphCorr(Ca,Pa,3,0,disRank);
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

% Set the powers of all local tests at rank 1 to 0 
power1(1,:)=0;power1(:,1)=0;power2(1,:)=0;power2(:,1)=0;power3(1,:)=0;power3(:,1)=0;
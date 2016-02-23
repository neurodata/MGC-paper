function [power1All, power2All, power3All, power4,power5,power6,power7]=CorrIndTestDim(type,n,dim,lim,rep1, rep2,noise,alpha,option)
% Author: Cencheng Shen
% Independence Tests for identifying dependency, with respect to increasing dimension at a fixed sample size.
% The output are the empirical powers of MGC{mcorr/dcorr/Mantel}, and global mcorr/dcorr/Mantel/HHG.
%
% Parameters:
% type specifies the type of distribution,
% n is the sample size, dim is the dimension,
% lim specifies the number of intervals in the sample size,
% rep1 specifies the number of MC-replicates for optimal scale estimation
% in MGC; if 0, the estimation step is skipped.
% rep2 specifies the number of MC-replicates for estimating the testing power,
% noise specifies the noise level, by default 1,
% alpha specifies the type 1 error level,
% option specifies whether each test statistic is calculated or not.
if nargin<5 || rep1<=0;
    rep1=1000;
end
if nargin<6 || rep2<=0;
    rep2=1000;
end
if nargin<7
    noise=0; % Default noise level
end
if nargin<8
    alpha=0.05; % Default type 1 error level
end
if nargin<9
    option=[1,1,1,1]; % Default option. Setting any to 0 to disable the calculation of MGC{mcorr/dcorr/Mantel} or HHG.
end

if lim==0
    dimRange=dim; % Test at dim only
    lim=1;
else
    dimRange=ceil(dim/lim):ceil(dim/lim):dim; % Test using dimension choices at the interval of ceil(dim/lim).
end
power1=zeros(1,lim);power2=zeros(1,lim);power3=zeros(1,lim);% Powers for MGC{mcorr/dcorr/Mantel}
power4=zeros(1,lim);power5=zeros(1,lim);power6=zeros(1,lim);% Powers for global mcorr/dcorr/Mantel.

% Run the independence test to first estimate the optimal scale of MGC
[~,~,~,~,neighborhoods]=IndependenceTestDim(type,n,dim,lim,rep1, noise,alpha); % Estimated optimal neighborhoods at each sample size. 

% Run the independence test again for the testing powers
[power1All, power2All, power3All, power7]=IndependenceTestDim(type,n,dim,lim,rep2, noise,alpha,option); % Powers for all local tests of mcorr/dcorr/Mantel, and HHG

% From the powers of all local tests, get the powers of MGC based on the optimal neighborhood estimation, and the powers of the respective global test
for i=1:lim
    tmp=power1All(:,:,i);
    power1(i)=tmp(neighborhoods(1,i));power4(i)=tmp(end,end);
    tmp=power2All(:,:,i);
    power2(i)=tmp(neighborhoods(2,i));power5(i)=tmp(end,end);
    tmp=power3All(:,:,i);
    power3(i)=tmp(neighborhoods(3,i));power6(i)=tmp(end,end);
end

% Save the results
pre1='../../Data/'; 
filename=strcat(pre1,'CorrIndTestDimType',num2str(type),'N',num2str(n),'Dim');
save(filename,'power1','power2','power3','power4','power5','power6','power7','type','n','rep1','rep2','lim','dim','noise','alpha','option','dimRange','neighborhoods','power1All','power2All','power3All');
% numRange=1:lim;
% plot(numRange,power1,'r.-',numRange,power2,'b.-',numRange,power3,'c.-',numRange,power4,'r.:',numRange,power5,'b.:',numRange,power6,'c.:',numRange,power7,'g.:','LineWidth',2);

function [power1, power2, power3, power4,neighbor]=IndependenceTestDim(type,n,dim,lim,rep, noise,alpha,option)
% This is an auxiliary function of the main function to calculate the powers of
% all local tests of mcorr/dcorr/Mantel, and the power of HHG.
%
% It first generates dependent and independent data from simulation
% distributions, then calculate the test statistics under the null and the
% alternative, and estimate the testing power of each method.
if nargin<8
    option=[1,1,1,0]; % Default option. Setting each entry to 0 to disable the calculation of local mcorr/dcorr/Mantel, or HHG.
end

if lim==0
    dimRange=dim; % Test at dim only
    lim=1;
else
    dimRange=ceil(dim/lim):ceil(dim/lim):dim; % Test using dimension choices at the interval of ceil(dim/lim).
end

% Store the test statistics under the null and the alternative
dCor1N=zeros(n,n,rep);dCor2N=zeros(n,n,rep);dCor3N=zeros(n,n,rep);
dCor1A=zeros(n,n,rep);dCor2A=zeros(n,n,rep);dCor3A=zeros(n,n,rep);
dCor4N=zeros(1,rep);dCor4A=zeros(1,rep);

% Powers
power1=zeros(n,n,lim);power2=zeros(n,n,lim);power3=zeros(n,n,lim);%Powers for all local tests of mcorr/dcorr/Mantel
power4=zeros(1,lim);% Powers for HHG
neighbor=zeros(3,lim); % Optimal neighborhoods for local mcorr/dcorr/Mantel

% Iterate through all dimension choices
for i=1:lim
    d=dimRange(i);
    % First generate null distribution, i.e., X and Y are independent
    % and estimate the distribution of the test statistics under the null
    for r=1:rep
        % Generate independent sample data
        [x,y]=CorrSampleGenerator(type,n,d,0, noise);
        % Form the distance matrix and calculate all test statistics under the null
        C=squareform(pdist(x));
        D=squareform(pdist(y));
        disRank=[disToRanks(C) disToRanks(D)];
        if option(1)~=0
            dCor1N(:,:,r)=LocalGraphCorr(C,D,1,disRank);
        end
        if option(2)~=0
            dCor2N(:,:,r)=LocalGraphCorr(C,D,2,disRank);
        end
        if option(3)~=0
            dCor3N(:,:,r)=LocalGraphCorr(C,D,3,disRank);
        end
        if option(4)~=0
            dCor4N(r)=HHG(C,D);
        end
    end
    
    % Then generate the alternative distribution, i.e., X and Y are dependent,
    % and estimate the distribution of the test statistics under the alternative
    for r=1:rep
        % Generate dependent sample data
        [x,y]=CorrSampleGenerator(type,n,d,1, noise);
        % Form the distance matrix and calculate all test statistics under the alternative
        C=squareform(pdist(x));
        D=squareform(pdist(y));
        disRank=[disToRanks(C) disToRanks(D)];
        if option(1)~=0
            dCor1A(:,:,r)=LocalGraphCorr(C,D,1,disRank);
        end
        if option(2)~=0
            dCor2A(:,:,r)=LocalGraphCorr(C,D,2,disRank);
        end
        if option(3)~=0
            dCor3A(:,:,r)=LocalGraphCorr(C,D,3,disRank);
        end
        if option(4)~=0
            dCor4A(r)=HHG(C,D);
        end
    end
    
    % Based on the emprical test statistics under the null and the alternative,
    % calculate the emprical powers and the best scale at type 1 error level alpha
    [power1(:,:,i),neighbor(1,i)]=calculatePower(dCor1N,dCor1A,alpha,rep);
    [power2(:,:,i),neighbor(2,i)]=calculatePower(dCor2N,dCor3A,alpha,rep);
    [power3(:,:,i),neighbor(3,i)]=calculatePower(dCor3N,dCor3A,alpha,rep);
    power4(i)=calculatePower(dCor4N,dCor4A,alpha,rep);
end

% Set the powers of all local tests at rank 0 to 0
% power1(1,:)=0;power1(:,1)=0;power2(1,:)=0;power2(:,1)=0;power3(1,:)=0;power3(:,1)=0;

function [power1,n1]=calculatePower(dCor1N,dCor1A,alpha,rep)
% An auxiliary function to estimate the power based on the distribution of
% the test statistic under the null and the alternative, and calculate the
% optimal scale from the estimated power
[n]=size(dCor1N,1);
power1=zeros(n,n);
if n>1
    for i=1:n;
        for j=1:n;
            dCorT=sort(dCor1N(i,j,:),'descend');
            cut1=dCorT(ceil(rep*alpha));
            power1(i,j)=mean(dCor1A(i,j,:)>cut1);
        end
    end
    n1=maxNeighbors(power1,dCor1N,dCor1A);
else
    dCorT=sort(dCor1N,'descend');
    cut1=dCorT(ceil(rep*alpha));
    power1=mean(dCor1A>cut1);
    n1=0;
end
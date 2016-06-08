function [power1, power2, power3, power4,power5,power6,power7]=CorrIndTest(type,n,dim,lim,rep1, rep2,noise,alpha,option)
% Author: Cencheng Shen
% Independence Tests for identifying dependency, with respect to increasing sample size at a fixed dimension.
% The output are the empirical powers of MGC{dcorr/mcorr/Mantel}, and global dcorr/mcorr/Mantel/HHG.
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
    noise=1; % Default noise level
end
if nargin<8
    alpha=0.05; % Default type 1 error level
end
if nargin<9
    option=[1,1,1,1]; % Default option. Setting any to 0 to disable the calculation of MGC{dcorr/mcorr/Mantel} or HHG.
end

if lim==0
    numRange=n; % Test at sample size n only
else
    numRange=ceil(n/lim):ceil(n/lim):n; % Test using sample sizes at the interval of ceil(n/lim).
end
lim=length(numRange);

power1=zeros(1,lim);power2=zeros(1,lim);power3=zeros(1,lim);% Powers for MGC{dcorr/mcorr/Mantel}
power4=zeros(1,lim);power5=zeros(1,lim);power6=zeros(1,lim);% Powers for global dcorr/mcorr/Mantel.

% Run the independence test to first estimate the optimal scale of MGC
[~,~,~,~,neighborhoods]=IndependenceTest(type,n,dim,lim,rep1, noise,alpha); % Estimated optimal neighborhoods at each sample size. 

% Run the independence test again for the testing powers
[power1All, power2All, power3All, power7]=IndependenceTest(type,n,dim,lim,rep2, noise,alpha,option); % Powers for all local tests of dcorr/mcorr/Mantel, and HHG

% From the powers of all local tests, get the powers of MGC based on the optimal neighborhood estimation, and the powers of the respective global test
for i=1:lim
    tmp=power1All(1:numRange(i),1:numRange(i),i);
    power1(i)=tmp(neighborhoods(1,i));
    tmp=tmp(tmp>0);
    power4(i)=tmp(end);
    tmp=power2All(1:numRange(i),1:numRange(i),i);
    power2(i)=tmp(neighborhoods(2,i));
    tmp=tmp(tmp>0);
    power5(i)=tmp(end);
    tmp=power3All(1:numRange(i),1:numRange(i),i);
    power3(i)=tmp(neighborhoods(3,i));
    tmp=tmp(tmp>0);
    power6(i)=tmp(end);
end

% Save the results
pre1='../../Data/Results/';
filename=strcat(pre1,'CorrIndTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
if type==0;
    filename=strcat(filename,'W',num2str(noise),'.mat');
end
save(filename,'power1','power2','power3','power4','power5','power6','power7','type','n','rep1','rep2','lim','dim','noise','alpha','option','numRange','neighborhoods','power1All','power2All','power3All');

function [power1, power2, power3, power4,neighbor]=IndependenceTest(type,n,dim,lim,rep, noise,alpha,option)
% This is an auxiliary function of the main function to calculate the powers of
% all local tests of dcorr/mcorr/Mantel, and the power of HHG.
%
% It first generates dependent and independent data from simulation
% distributions, then calculate the test statistics under the null and the
% alternative, and estimate the testing power of each method.
if nargin<8
    option=[1,1,1,0]; % Default option. Setting each entry to 0 to disable the calculation of local dcorr/mcorr/Mantel, or HHG.
end

if lim==0
    numRange=n; % Test at sample size n only
else
    numRange=ceil(n/lim):ceil(n/lim):n; % Test using sample sizes at the interval of ceil(n/lim).
end
lim=length(numRange);
d=dim;

% Store the test statistics under the null and the alternative
dCor1N=zeros(n,n,rep);dCor2N=zeros(n,n,rep);dCor3N=zeros(n,n,rep);
dCor1A=zeros(n,n,rep);dCor2A=zeros(n,n,rep);dCor3A=zeros(n,n,rep);
dCor4N=zeros(1,rep);dCor4A=zeros(1,rep);
% Store the dependent and independent data
DataN=zeros(n,2*n,rep);DataA=zeros(n,2*n,rep);

% Powers
power1=zeros(n,n,lim);power2=zeros(n,n,lim);power3=zeros(n,n,lim);% Powers for all local tests of dcorr/mcorr/Mantel
power4=zeros(1,lim);% Powers for HHG
neighbor=zeros(3,lim); % Optimal neighborhoods for local dcorr/mcorr/Mantel

for r=1:rep
    % Generate independent sample data and form the distance matrices
    [x,y]=CorrSampleGenerator(type,n,d,0, noise);
    C=squareform(pdist(x));
    D=squareform(pdist(y));
    DataN(:,:,r)=[C D];
    
    % Generate dependent sample data and form the distance matrices
    [x,y]=CorrSampleGenerator(type,n,d,1, noise);
    C=squareform(pdist(x));
    D=squareform(pdist(y));
    DataA(:,:,r)=[C D];
end

% Iterate through all sample sizes
for i=1:lim
    nn=numRange(i);
    % First estimate the distribution of the test statistics under the null
    for r=1:rep
        C=DataN(1:nn,1:nn,r);
        D=DataN(1:nn,n+1:n+nn,r);
        if option(1)~=0
            tmp=LocalCorr(C,D,1);
            dCor1N(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(2)~=0
            tmp=LocalCorr(C,D,2);
            dCor2N(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(3)~=0
            tmp=LocalCorr(C,D,3);
            dCor3N(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(4)~=0
            dCor4N(r)=HHG(C,D);
        end
    end
    
    % Then estimate the distribution of the test statistics under the alternative
    for r=1:rep
        C=DataA(1:nn,1:nn,r);
        D=DataA(1:nn,n+1:n+nn,r);
        if option(1)~=0
            tmp=LocalCorr(C,D,1);
            dCor1A(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(2)~=0
            tmp=LocalCorr(C,D,2);
            dCor2A(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(3)~=0
            tmp=LocalCorr(C,D,3);
            dCor3A(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(4)~=0
            dCor4A(r)=HHG(C,D);
        end
    end
    
    % Based on the emprical test statistics under the null and the alternative,
    % calculate the emprical powers and the best scale at type 1 error level alpha
    [power1(1:nn,1:nn,i),neighbor(1,i)]=calculatePower(dCor1N(1:nn,1:nn,:),dCor1A(1:nn,1:nn,:),alpha,rep);
    [power2(1:nn,1:nn,i),neighbor(2,i)]=calculatePower(dCor2N(1:nn,1:nn,:),dCor2A(1:nn,1:nn,:),alpha,rep);
    [power3(1:nn,1:nn,i),neighbor(3,i)]=calculatePower(dCor3N(1:nn,1:nn,:),dCor3A(1:nn,1:nn,:),alpha,rep);
    power4(i)=calculatePower(dCor4N,dCor4A,alpha,rep);
end

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
    power1(1,:)=0;power1(:,1)=0; % Set the powers of all local tests at rank 0 to 0
    n1=maxNeighbors(power1,0,dCor1N,dCor1A);
else
    dCorT=sort(dCor1N,'descend');
    cut1=dCorT(ceil(rep*alpha));
    power1=mean(dCor1A>cut1);
    n1=0;
end
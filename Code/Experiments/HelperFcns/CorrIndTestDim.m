function [powerMGC,powerDLocal,powerMLocal,powerPLocal,powerD,powerM,powerP,powerHHG]=CorrIndTestDim(type,n,dim,lim,rep1, rep2,noise,alpha,option)
% Author: Cencheng Shen
% Independence Tests for identifying dependency, with respect to increasing dimension at a fixed sample size.
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
    noise=0; % Default noise level
end
if nargin<8
    alpha=0.05; % Default type 1 error level
end
if nargin<9
    option=[1,1,1,1]; % Default option. Setting any to 0 to disable the calculation of MGC{dcorr/mcorr/Mantel} or HHG.
end

if lim==0
    dimRange=dim; % Test at dim only
    lim=1;
else
    dimRange=ceil(dim/lim):ceil(dim/lim):dim; % Test using dimension choices at the interval of ceil(dim/lim).
end
if ceil(dim/lim)~=1
    dimRange=[1 dimRange];
    lim=length(dimRange);
end
powerMGCD=zeros(1,lim);powerMGCM=zeros(1,lim);powerMGCP=zeros(1,lim);% Powers for MGC{dcorr/mcorr/Mantel}
powerD=zeros(1,lim);powerM=zeros(1,lim);powerP=zeros(1,lim);% Powers for global dcorr/mcorr/Mantel.

% Run the independence test to first estimate the optimal scale of MGC
[~,~,~,~,~,neighborhoods]=IndependenceTestDim(type,n,dimRange,lim,rep1, noise,alpha); % Estimated optimal neighborhoods at each sample size.

% Run the independence test again for the testing powers
[powerMGC,powerDLocal, powerMLocal, powerPLocal, powerHHG]=IndependenceTestDim(type,n,dimRange,lim,rep2, noise,alpha,option); % Powers for all local tests of dcorr/mcorr/Mantel, and HHG

% From the powers of all local tests, get the powers of MGC based on the optimal neighborhood estimation, and the powers of the respective global test
for i=1:lim
    if option(1)~=0
        tmp=powerDLocal(:,:,i);
        powerMGCD(i)=tmp(neighborhoods(1,i));
        tmp=tmp(tmp>0);
        powerD(i)=tmp(end);
    end
    if option(2)~=0
        tmp=powerMLocal(:,:,i);
        powerMGCM(i)=tmp(neighborhoods(2,i));
        tmp=tmp(tmp>0);
        powerM(i)=tmp(end);
    end
    if option(3)~=0
        tmp=powerPLocal(:,:,i);
        powerMGCP(i)=tmp(neighborhoods(3,i));
        tmp=tmp(tmp>0);
        powerP(i)=tmp(end);
    end
end

% Save the results
%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-3));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

pre1=strcat(rootDir,'Data/Results/');% The folder to save figures
filename=strcat(pre1,'CorrIndTestDimType',num2str(type),'N',num2str(n),'Dim');
save(filename,'powerMGCD','powerMGCM','powerMGCP','powerD','powerM','powerP','powerHHG','powerMGC','type','n','rep1','rep2','lim','dim','noise','alpha','option','dimRange','neighborhoods','powerDLocal','powerMLocal','powerPLocal');

function [powerMGC,powerD, powerM, powerP, powerHHG,neighbor]=IndependenceTestDim(type,n,dimRange,lim,rep, noise,alpha,option)
% This is an auxiliary function of the main function to calculate the powers of
% all local tests of dcorr/mcorr/Mantel, and the power of HHG.
%
% It first generates dependent and independent data from simulation
% distributions, then calculate the test statistics under the null and the
% alternative, and estimate the testing power of each method.
if nargin<8
    option=[1,1,1,0]; % Default option. Setting each entry to 0 to disable the calculation of local dcorr/mcorr/Mantel, or HHG.
end

% Store the test statistics under the null and the alternative
dCorDN=zeros(n,n,rep);dCorMN=zeros(n,n,rep);dCorPN=zeros(n,n,rep);
dCorDA=zeros(n,n,rep);dCorMA=zeros(n,n,rep);dCorPA=zeros(n,n,rep);
dCorHHGN=zeros(1,rep);dCorHHGA=zeros(1,rep);
dCorMGCN=zeros(1,rep);dCorMGCA=zeros(1,rep);

% Powers
powerD=zeros(n,n,lim);powerM=zeros(n,n,lim);powerP=zeros(n,n,lim);% Powers for all local tests of dcorr/mcorr/Mantel
powerHHG=zeros(1,lim);% Powers for HHG
powerMGC=zeros(1,lim);% Powers for HHG
neighbor=zeros(3,lim); % Optimal neighborhoods for local dcorr/mcorr/Mantel

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
        if option(1)~=0
            tmp=LocalCorr(C,D,'dcor');
            dCorDN(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(2)~=0
            tmp=LocalCorr(C,D,'mcor');
            dCorMN(1:size(tmp,1),1:size(tmp,2),r)=tmp;
            dCorMGCN(r)=MGCSampleStat(tmp);
        end
        if option(3)~=0
            tmp=LocalCorr(C,D,'mantel');
            dCorPN(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(4)~=0
            dCorHHGN(r)=HHG(C,D);
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
        if option(1)~=0
            tmp=LocalCorr(C,D,'dcor');
            dCorDA(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(2)~=0
            tmp=LocalCorr(C,D,'mcor');
            dCorMA(1:size(tmp,1),1:size(tmp,2),r)=tmp;
            dCorMGCA(r)=MGCSampleStat(tmp);
        end
        if option(3)~=0
            tmp=LocalCorr(C,D,'mantel');
            dCorPA(1:size(tmp,1),1:size(tmp,2),r)=tmp;
        end
        if option(4)~=0
            dCorHHGA(r)=HHG(C,D);
        end
    end
    
    % Based on the emprical test statistics under the null and the alternative,
    % calculate the emprical powers and the best scale at type 1 error level alpha
    [powerD(:,:,i),neighbor(1,i)]=calculatePower(dCorDN,dCorDA,alpha,rep);
    [powerM(:,:,i),neighbor(2,i)]=calculatePower(dCorMN,dCorMA,alpha,rep);
    [powerP(:,:,i),neighbor(3,i)]=calculatePower(dCorPN,dCorPA,alpha,rep);
    powerHHG(i)=calculatePower(dCorHHGN,dCorHHGA,alpha,rep);
    powerMGC(i)=calculatePower(dCorMGCN,dCorMGCA,alpha,rep);
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
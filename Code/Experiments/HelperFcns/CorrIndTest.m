function [powerMGC,powerDLocal,powerMLocal,powerPLocal,powerD,powerM,powerP,powerHHG,powerHSIC]=CorrIndTest(type,n,dim,lim,rep,noise,alpha,option)
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
if nargin<5 || rep<=0;
    rep=1000;
end
if nargin<6
    noise=1; % Default noise level
end
if nargin<7
    alpha=0.05; % Default type 1 error level
end
if nargin<8
    option=[1,1,1,1,1,1,1,1,1]; % Default option. Setting any to 0 to disable the calculation of MGC{mcorr/dcorr/Mantel} or HHG.
end

if lim==0
    numRange=n; % Test at sample size n only
else
    numRange=ceil(n/lim):ceil(n/lim):n; % Test using sample sizes at the interval of ceil(n/lim).
end
lim=length(numRange);

powerMGCD=zeros(1,lim);powerMGCM=zeros(1,lim);powerMGCP=zeros(1,lim);% Powers for MGC{dcorr/mcorr/Mantel}
powerD=zeros(1,lim);powerM=zeros(1,lim);powerP=zeros(1,lim);% Powers for global dcorr/mcorr/Mantel.
powerHSIC=zeros(1,lim);powerHHG=zeros(1,lim);
powerCorr=zeros(1,lim);powerSpearman=zeros(1,lim);powerKendall=zeros(1,lim);powerMIC=zeros(1,lim);
neighborhoods=ones(3,lim);

% Run the independence test to first estimate the optimal scale of MGC
% [~,~,~,~,~,~,neighborhoods]=IndependenceTest(type,numRange,dim,lim,rep1, noise,alpha); % Estimated optimal neighborhoods at each sample size.

% Run the independence test again for the testing powers
[powerMGC,powerMGC2, powerMGC3,powerDLocal, powerMLocal, powerPLocal, powerHHG,powerHSIC,powerCorr,powerSpearman,powerKendall,powerMIC]=IndependenceTest(type,numRange,dim,lim,rep, noise,alpha,option); % Powers for all local tests of dcorr/mcorr/Mantel, and HHG

% From the powers of all local tests, get the powers of MGC based on the optimal neighborhood estimation, and the powers of the respective global test
for i=1:lim
    if option(1)~=0
        tmp=powerMLocal(1:numRange(i),1:numRange(i),i);
        powerMGCM(i)=max(max(tmp));
%         powerMGCM(i)=tmp(neighborhoods(2,i));
%         tmp=tmp(tmp>0);
        powerM(i)=tmp(end);
    end
    if option(2)~=0
        tmp=powerDLocal(1:numRange(i),1:numRange(i),i);
        powerMGCD(i)=max(max(tmp));
%         powerMGCD(i)=tmp(neighborhoods(1,i));
%         tmp=tmp(tmp>0);
        powerD(i)=tmp(end);
    end
    if option(3)~=0
        tmp=powerPLocal(1:numRange(i),1:numRange(i),i);
        powerMGCP(i)=max(max(tmp));
%         powerMGCP(i)=tmp(neighborhoods(3,i));
%         tmp=tmp(tmp>0);
        powerP(i)=tmp(end);
    end
end

% Save the results
%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-3));
addpath(genpath(strcat(rootDir,'Code/')));

pre1=strcat(rootDir,'Data/Results/');% The folder to save figures
filename=strcat(pre1,'CorrIndTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
if type==0;
    filename=strcat(filename,'W',num2str(noise),'.mat');
end
save(filename,'powerCorr','powerMGC2', 'powerMGC3','powerSpearman','powerKendall','powerMIC','powerMGCD','powerMGCM','powerMGCP','powerD','powerM','powerP','powerHHG','powerHSIC','powerMGC','type','n','rep','lim','dim','noise','alpha','option','numRange','neighborhoods','powerDLocal','powerMLocal','powerPLocal');

function [powerMGC,powerMGC2, powerMGC3,powerD, powerM, powerP, powerHHG,powerHSIC,powerCorr,powerSpearman,powerKendall,powerMIC,neighbor]=IndependenceTest(type,numRange,dim,lim,rep, noise,alpha,option)
% This is an auxiliary function of the main function to calculate the powers of
% all local tests of dcorr/mcorr/Mantel, and the power of HHG.
%
% It first generates dependent and independent data from simulation
% distributions, then calculate the test statistics under the null and the
% alternative, and estimate the testing power of each method.
if nargin<8
    option=[1,1,1,1,1,1,1,1,1]; % Default option. Setting each entry to 0 to disable the calculation of local dcorr/mcorr/Mantel, or HHG.
end

d=dim;
n=numRange(end);

% Store the test statistics under the null and the alternative
dCorDN=zeros(n,n,rep);dCorMN=zeros(n,n,rep);dCorPN=zeros(n,n,rep);
dCorDA=zeros(n,n,rep);dCorMA=zeros(n,n,rep);dCorPA=zeros(n,n,rep);
dCorHHGN=zeros(1,rep);dCorHHGA=zeros(1,rep);
dCorHSICN=zeros(1,rep);dCorHSICA=zeros(1,rep);
dCorMGCN=zeros(1,rep);dCorMGCA=zeros(1,rep);
dCorMGC2N=zeros(1,rep);dCorMGC2A=zeros(1,rep);
dCorMGC3N=zeros(1,rep);dCorMGC3A=zeros(1,rep);
dCorCorrN=zeros(1,rep);dCorCorrA=zeros(1,rep);
dCorSpearmanN=zeros(1,rep);dCorSpearmanA=zeros(1,rep);
dCorKendallN=zeros(1,rep);dCorKendallA=zeros(1,rep);
dCorMICN=zeros(1,rep);dCorMICA=zeros(1,rep);
% Store the dependent and independent data
DataN=zeros(d,2*n,rep);DataA=zeros(d,2*n,rep);

% Powers
powerD=zeros(n,n,lim);powerM=zeros(n,n,lim);powerP=zeros(n,n,lim);% Powers for all local tests of dcorr/mcorr/Mantel
powerHHG=zeros(1,lim);% Powers for HHG
powerHSIC=zeros(1,lim);% Powers for HHG
powerMGC=zeros(1,lim);% Powers for HHG
powerMGC2=zeros(1,lim);% Powers for HHG
powerMGC3=zeros(1,lim);% Powers for HHG
powerCorr=zeros(1,lim);% Powers for HHG
powerSpearman=zeros(1,lim);% Powers for HHG
powerKendall=zeros(1,lim);% Powers for HHG
powerMIC=zeros(1,lim);% Powers for HHG
neighbor=zeros(3,lim); % Optimal neighborhoods for local dcorr/mcorr/Mantel

for r=1:rep
    % Generate independent sample data and form the distance matrices
    [x,y]=CorrSampleGenerator(type,n,d,0, noise);
%     C=squareform(pdist(x));
%     D=squareform(pdist(y));
    DataN(:,:,r)=[x' y'];
    
    % Generate dependent sample data and form the distance matrices
    [x,y]=CorrSampleGenerator(type,n,d,1, noise);
%     C=squareform(pdist(x));
%     D=squareform(pdist(y));
    DataA(:,:,r)=[x' y'];
end

% Iterate through all sample sizes
for i=1:lim
    nn=numRange(i);
    % First estimate the distribution of the test statistics under the null
    for r=1:rep
        x=DataN(:,1:nn,r)';
        y=DataN(:,n+1:n+nn,r)';
        C=squareform(pdist(x));
        D=squareform(pdist(y));
        if option(1)~=0
            [stat, tmp, ~]=MGCSampleStat(C,D,'mcor');     
            dCorMN(1:size(tmp,1),1:size(tmp,2),r)=tmp;
            dCorMGCN(r)=stat;
        end
        if option(2)~=0
            [stat, tmp, ~]=MGCSampleStat(C,D,'dcor');     
            dCorDN(1:size(tmp,1),1:size(tmp,2),r)=tmp;
            dCorMGC2N(r)=stat;
        end
        if option(3)~=0
            [stat, tmp, ~]=MGCSampleStat(C,D,'mantel');     
            dCorPN(1:size(tmp,1),1:size(tmp,2),r)=tmp;
            dCorMGC3N(r)=stat;
        end
        if option(4)~=0
            dCorHHGN(r)=HHG(C,D);
        end
        if option(5)~=0
            dCorHSICN(r)=HSIC(x,y);
        end
        if option(6)~=0
            dCorCorrN(r)=corr(x,y);
        end
        if option(7)~=0
            dCorSpearmanN(r)=corr(x,y,'type','Spearman');
        end
        if option(8)~=0
            dCorKendallN(r)=corr(x,y,'type','Kendall');
        end
        if option(9)~=0
            an=mine(x',y');
            dCorMICN(r)=an.mic;
        end
    end
    
    % Then estimate the distribution of the test statistics under the alternative
    for r=1:rep
        x=DataA(:,1:nn,r)';
        y=DataA(:,n+1:n+nn,r)';
        C=squareform(pdist(x));
        D=squareform(pdist(y));
        if option(1)~=0
            [stat, tmp, ~]=MGCSampleStat(C,D,'mcor');     
            dCorMA(1:size(tmp,1),1:size(tmp,2),r)=tmp;
            dCorMGCA(r)=stat;
        end
        if option(2)~=0
            [stat, tmp, ~]=MGCSampleStat(C,D,'dcor');     
            dCorDA(1:size(tmp,1),1:size(tmp,2),r)=tmp;
            dCorMGC2A(r)=stat;
        end
        if option(3)~=0
            [stat, tmp, ~]=MGCSampleStat(C,D,'mantel');     
            dCorPA(1:size(tmp,1),1:size(tmp,2),r)=tmp;
            dCorMGC3A(r)=stat;
        end
        if option(4)~=0
            dCorHHGA(r)=HHG(C,D);
        end
        if option(5)~=0
            dCorHSICA(r)=HSIC(x,y);
        end
        if option(6)~=0
            dCorCorrA(r)=corr(x,y);
        end
        if option(7)~=0
            dCorSpearmanA(r)=corr(x,y,'type','Spearman');
        end
        if option(8)~=0
            dCorKendallA(r)=corr(x,y,'type','Kendall');
        end
        if option(9)~=0
            an=mine(x',y');
            dCorMICA(r)=an.mic;
        end
    end
    
    % Based on the emprical test statistics under the null and the alternative,
    % calculate the emprical powers and the best scale at type 1 error level alpha
    [powerD(1:nn,1:nn,i),neighbor(1,i)]=calculatePower(dCorDN(1:nn,1:nn,:),dCorDA(1:nn,1:nn,:),alpha,rep);
    [powerM(1:nn,1:nn,i),neighbor(2,i)]=calculatePower(dCorMN(1:nn,1:nn,:),dCorMA(1:nn,1:nn,:),alpha,rep);
    [powerP(1:nn,1:nn,i),neighbor(3,i)]=calculatePower(dCorPN(1:nn,1:nn,:),dCorPA(1:nn,1:nn,:),alpha,rep);
    powerHHG(i)=calculatePower(dCorHHGN,dCorHHGA,alpha,rep);
    powerHSIC(i)=calculatePower(dCorHSICN,dCorHSICA,alpha,rep);
    powerMGC(i)=calculatePower(dCorMGCN,dCorMGCA,alpha,rep);
    powerMGC2(i)=calculatePower(dCorMGC2N,dCorMGC2A,alpha,rep);
    powerMGC3(i)=calculatePower(dCorMGC3N,dCorMGC3A,alpha,rep);
    powerCorr(i)=calculatePower(dCorCorrN,dCorCorrA,alpha,rep);
    powerSpearman(i)=calculatePower(dCorSpearmanN,dCorSpearmanA,alpha,rep);
    powerKendall(i)=calculatePower(dCorKendallN,dCorKendallA,alpha,rep);
    powerMIC(i)=calculatePower(dCorMICN,dCorMICA,alpha,rep);
%     figure
%     hold on
%     [f,x]=ksdensity(reshape(dCorMN(nn,nn,:),1,rep));
%     plot(x,f,'r+')
%     [f,x]=ksdensity(reshape(dCorMA(nn,nn,:),1,rep));
%     plot(x,f,'r.')
%     [f,x]=ksdensity(dCorMGCN);
%     plot(x,f,'g-')
%     [f,x]=ksdensity(dCorMGCA);
%     plot(x,f,'g-')
%     s=strcat(num2str(powerMGC(i)),',',num2str(powerM(nn,nn,i)),',',num2str(max(max(powerM(:,:,i)))));
%     title(s);
%     hold off
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
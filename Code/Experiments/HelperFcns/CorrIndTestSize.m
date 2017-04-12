function [size]=CorrIndTestSize(type,nMax,dim,thres,rep,noise,alpha,option)
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
if nargin<4
    thres=0.85;
end
if nargin<5 || rep<=0
    rep=1000;
end
if nargin<6
    noise=1; % Default noise level
end
if nargin<7
    alpha=0.05; % Default type 1 error level
end
if nargin<8
    option=0; % Default option, 0 for oracle MGC, 1 for sample mgc, 2 for dcorr, 3 for mcorr, 4 for mantel, 5 for HHG
end

nRange=10:10:nMax;
ln=length(nRange);
% Run the independence test to first estimate the optimal scale of MGC
% Run the independence test again for the testing powers
ind=min(ln,5);
if type==12 || type ==13
    ind=ln;
end
indMin=ind;
power=0;
repp=0;
mm=abs(power-thres);
while mm>0.01 && repp<10
    repp=repp+1;
    n=nRange(ind);
    [power]=IndependenceTest(type,n,dim,rep, noise,alpha,option); % Powers for all local tests of dcorr/mcorr/Mantel, and HHG
    if option==0
        [~,neighborhoods]=IndependenceTest(type,n,dim,rep, noise,alpha,option); % Estimated optimal neighborhoods at each sample size.
        power=power(neighborhoods);
    end
    tmp=abs(power-thres);
    if tmp<mm
        mm=tmp;
        indMin=ind;
    end
    if power<thres && ind==ln
        indMin=ind;
        break;
    end
    if power>thres && ind==1
        indMin=ind;
        break;
    end
    if power<0.1
        ind=min(4*ind,ln);
        continue;
    end
    if power>0.99
        ind=max(1,ceil(ind/2));
        continue;
    end
    ind=min(ind+ceil((thres-power)/0.1),ln);
end

size=nRange(indMin);

% % Save the results
% %%% File path searching
% fpath = mfilename('fullpath');
% fpath=strrep(fpath,'\','/');
% findex=strfind(fpath,'/');
% rootDir=fpath(1:findex(end-3));
% addpath(genpath(strcat(rootDir,'Code/')));
% 
% pre1=strcat(rootDir,'Data/Results/');% The folder to save figures
% filename=strcat(pre1,'CorrIndTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
% if type==0
%     filename=strcat(filename,'W',num2str(noise),'.mat');
% end
% save(filename,'powerMGCD','powerMGCM','powerMGCP','powerD','powerM','powerP','powerHHG','powerMGC','type','n','rep1','rep2','lim','dim','noise','alpha','option','numRange','neighborhoods','powerDLocal','powerMLocal','powerPLocal');

function [power,scale]=IndependenceTest(type,n,d,rep, noise,alpha,option)
% This is an auxiliary function of the main function to calculate the powers of
% all local tests of dcorr/mcorr/Mantel, and the power of HHG.
%
% It first generates dependent and independent data from simulation
% distributions, then calculate the test statistics under the null and the
% alternative, and estimate the testing power of each method.

% Store the test statistics under the null and the alternative
if option==0
    dCorN=zeros(n,n,rep);dCorA=zeros(n,n,rep); 
else
    dCorN=zeros(1,rep);dCorA=zeros(1,rep);
end
% Powers
switch option
    case {0,1,2}
        opt='mcor';
    case 3
        opt='dcor';
    case 4
        opt='mantel';
end
for r=1:rep
    % Generate independent sample data and form the distance matrices
    [x,y]=CorrSampleGenerator(type,n,d,0, noise);
    C=squareform(pdist(x));
    D=squareform(pdist(y));
    if option==5
        dCorN(r)=HHG(C,D);
    else
        tmp=MGCLocalCorr(C,D,opt);
        if option==1
            tmp=MGCSampleStat(tmp);
        end
        if option~=0
            dCorN(r)=tmp(end);
        else
            dCorN(:,:,r)=tmp;
        end
    end
    
    % Generate dependent sample data and form the distance matrices
    [x,y]=CorrSampleGenerator(type,n,d,1, noise);
    C=squareform(pdist(x));
    D=squareform(pdist(y));
    if option==5
        dCorA(r)=HHG(C,D);
    else
        tmp=MGCLocalCorr(C,D,opt);
        if option==1
            tmp=MGCSampleStat(tmp);
        end
        if option~=0
            dCorA(r)=tmp(end);
        else
            dCorA(:,:,r)=tmp;
        end
    end
end

% Based on the emprical test statistics under the null and the alternative,
% calculate the emprical powers and the best scale at type 1 error level alpha
[power,scale]=calculatePower(dCorN,dCorA,alpha,rep);

function [power1,n1]=calculatePower(dCor1N,dCor1A,alpha,rep)
% An auxiliary function to estimate the power based on the distribution of
% the test statistic under the null and the alternative, and calculate the
% optimal scale from the estimated power
[n]=size(dCor1N,1);
power1=zeros(n,n);
if n>1
    for i=1:n
        for j=1:n
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
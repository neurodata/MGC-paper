function [power1, power2, power3, power4,power5,power6,power7]=CorrIndTest(type,n,dim,lim,rep1, rep2,noise,option)
% Author: Cencheng Shen
% Independence Tests for identifying dependency, with respect to increasing sample size at a fixed dimension.
% The output are the empirical powers of MGC by mcorr/dcorr/Mantel, and global mcorr/dcorr/Mantel/HHG.
%
% Parameters:
% type specifies the type of distribution,
% n is the sample size, dim is the dimension,
% lim specifies the number of intervals in the sample size,
% rep1 specifies the number of MC-replicates for optimal scale estimation
% in MGC; if 0, the estimation step is skipped and the powers of all local
% tests are returned instead.
% rep2 specifies the number of MC-replicates for estimating the testing power,
% noise specifies the noise level, by default 1,
% option specifies whether each test statistic is calculated or not
if nargin<7
    noise=1; % Default noise level
end
if nargin<8
    option=[1,1,1,1,1,1,1]; % Default option. Setting any to 0 to disable the calculation of MGC by mcorr/dcorr/Mantel, global mcorr/dcorr/Mantel, HHG, in order.
end

if lim==0
    numRange=n; % Test at sample size n only
else
    numRange=ceil(n/lim):ceil(n/lim):n; % Test using sample sizes at the interval of ceil(n/lim).
end
lim=length(numRange);

power1=zeros(1,lim);power2=zeros(1,lim);power3=zeros(1,lim);% Powers for MGC by mcorr/dcorr/Mantel
power4=zeros(1,lim);power5=zeros(1,lim);power6=zeros(1,lim);power7=zeros(1,lim);% Powers for mcorr/dcorr/Mantel/HHG.
neighborhoods=zeros(3,lim); % Estimated optimal neighborhoods at each sample size. At 0, MGC of all scales are calculated

% Run the independence test to first estimate the optimal scale of MGC
if rep1~=0
    [p1,p2,p3]=IndependenceTest(type,n,dim,lim,rep1, noise);
    % Find the best scale, and take the last one when there exists more than one optimal scale
    for i=1:lim
        [neighbor1]=verifyNeighbors(1-p1(1:numRange(i),1:numRange(i),i));
        neighborhoods(1,i)=neighbor1(end);
        
        [neighbor2]=verifyNeighbors(1-p2(1:numRange(i),1:numRange(i),i));
        neighborhoods(2,i)=neighbor2(end);
        
        [neighbor3]=verifyNeighbors(1-p3(1:numRange(i),1:numRange(i),i));
        neighborhoods(3,i)=neighbor3(end);
    end
end

% Run the independence test again for the testing powers
[power1, power2, power3, power4,power5,power6,power7]=IndependenceTest(type,n,dim,lim,rep2, noise,option,neighborhoods);
% Save the results
if rep1==0
    tmpC='All';
else
    tmpC='';
end
filename=strcat('CorrIndTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim),tmpC);
save(filename,'power1','power2','power3','power4','power5','power6','power7','type','n','rep1','rep2','lim','dim','noise','option','numRange','neighborhoods');

%Plot
% numRange=1:lim;
% plot(numRange,power1,'ro-',numRange,power2,'bx-',numRange,power3,'c+-',numRange,power4,'r.:',numRange,power5,'b.:',numRange,power6,'c.:',numRange,power7,'g.:','LineWidth',2);
%

function [power1, power2, power3, power4,power5,power6,power7]=IndependenceTest(type,n,dim,lim,rep, noise,option,neighborhoods)
% Author: Cencheng Shen
% This is an auxiliary function of the main test.
%
% It first generates dependent and independent data from simulation
% distributions, calculate the test statistics under the null and the
% alternative, then estimate the testing power of each method.
if nargin<6
    noise=1; % Default noise level
end
if nargin<7
    option=[1,1,1,0,0,0,0]; % Default option to only consider MGC.
end
if nargin<8 || size(neighborhoods,1)<3
    neighborhoods=zeros(3,lim); % Default neighborhood to computes MGC of all scales
end
alpha=0.05; % Type 1 error level

if lim==0
    numRange=n; % Test at sample size n only
else
    numRange=ceil(n/lim):ceil(n/lim):n; % Test using sample sizes at the interval of ceil(n/lim).
end
lim=length(numRange);
d=dim;

% Store the test statistics under the null and the alternative
if norm(neighborhoods,'fro')==0
    indD=[n,n];
else
    indD=1;
end
dCor1N=zeros([indD,rep]);dCor2N=zeros([indD,rep]);dCor3N=zeros([indD,rep]);
dCor1A=zeros([indD,rep]);dCor2A=zeros([indD,rep]);dCor3A=zeros([indD,rep]);
dCor4N=zeros(1,rep);dCor4A=zeros(1,rep);dCor5N=zeros(1,rep);dCor5A=zeros(1,rep);
dCor6N=zeros(1,rep);dCor6A=zeros(1,rep);dCor7N=zeros(1,rep);dCor7A=zeros(1,rep);
% Store the dependent and independent data
DataN=zeros(n,2*n,rep);DataA=zeros(n,2*n,rep);

% Powers
power1=zeros([indD,lim]);power2=zeros([indD,lim]);power3=zeros([indD,lim]);% Powers for MGC by mcorr/dcorr/Mantel
power4=zeros(1,lim);power5=zeros(1,lim);power6=zeros(1,lim);power7=zeros(1,lim);% Powers for mcorr/dcorr/Mantel/HHG

for r=1:rep
    % Generate independent sample data and form the distance matrices
    [x,y]=CorrSampleGenerator(type,n,d,0, noise);
    C1=squareform(pdist(x));
    P1=squareform(pdist(y));
    DataN(:,:,r)=[C1 P1];
    
    % Generate dependent sample data and form the distance matrices
    [x,y]=CorrSampleGenerator(type,n,d,1, noise);
    C1=squareform(pdist(x));
    P1=squareform(pdist(y));
    DataA(:,:,r)=[C1 P1];
end

% Iterate through all sample sizes
for i=1:lim
    nn=numRange(i);
    % First estimate the distribution of the test statistics under the null
    for r=1:rep
        C=DataN(1:nn,1:nn,r);
        P=DataN(1:nn,n+1:n+nn,r);
        disRank=[disToRanks(C) disToRanks(P)];
        if option(1)~=0
            tmp=LocalGraphCorr(C,P,1,neighborhoods(1,i),disRank);
            if neighborhoods(1,i)==0;
                dCor1N(1:nn,1:nn,r)=tmp;
            else
                dCor1N(r)=tmp;
            end
        end
        if option(2)~=0
            tmp=LocalGraphCorr(C,P,2,neighborhoods(2,i),disRank);
            if neighborhoods(2,i)==0;
                dCor2N(1:nn,1:nn,r)=tmp;
            else
                dCor2N(r)=tmp;
            end
        end
        if option(3)~=0
            tmp=LocalGraphCorr(C,P,3,neighborhoods(3,i),disRank);
            if neighborhoods(3,i)==0;
                dCor3N(1:nn,1:nn,r)=tmp;
            else
                dCor3N(r)=tmp;
            end
        end
        if option(4)~=0
            dCor4N(r)=LocalGraphCorr(C,P,1,nn^2,disRank);
        end
        if option(5)~=0
            dCor5N(r)=LocalGraphCorr(C,P,2,nn^2,disRank);
        end
        if option(6)~=0
            dCor6N(r)=LocalGraphCorr(C,P,3,nn^2,disRank);
        end
        if option(7)~=0
            dCor7N(r)=HHG(C,P);
        end
    end
    
    % Then estimate the distribution of the test statistics under the alternative
    for r=1:rep
        C=DataA(1:nn,1:nn,r);
        P=DataA(1:nn,n+1:n+nn,r);
        disRank=[disToRanks(C) disToRanks(P)];
        if option(1)~=0
            tmp=LocalGraphCorr(C,P,1,neighborhoods(1,i),disRank);
            if neighborhoods(1,i)==0;
                dCor1A(1:nn,1:nn,r)=tmp;
            else
                dCor1A(r)=tmp;
            end
        end
        if option(2)~=0
            tmp=LocalGraphCorr(C,P,2,neighborhoods(2,i),disRank);
            if neighborhoods(2,i)==0;
                dCor2A(1:nn,1:nn,r)=tmp;
            else
                dCor2A(r)=tmp;
            end
        end
        if option(3)~=0
            tmp=LocalGraphCorr(C,P,3,neighborhoods(3,i),disRank);
            if neighborhoods(3,i)==0;
                dCor3A(1:nn,1:nn,r)=tmp;
            else
                dCor3A(r)=tmp;
            end
        end
        if option(4)~=0
            dCor4A(r)=LocalGraphCorr(C,P,1,nn^2,disRank);
        end
        if option(5)~=0
            dCor5A(r)=LocalGraphCorr(C,P,2,nn^2,disRank);
        end
        if option(6)~=0
            dCor6A(r)=LocalGraphCorr(C,P,3,nn^2,disRank);
        end
        if option(7)~=0
            dCor7A(r)=HHG(C,P);
        end
    end
    
    % Based on the emprical test statistics under the null and the alternative,
    % calculate the emprical powers at type 1 error level alpha=0.05
    power1=calculatePower(power1,i,numRange,dCor1N,dCor1A,alpha,rep);
    power2=calculatePower(power2,i,numRange,dCor2N,dCor2A,alpha,rep);
    power3=calculatePower(power3,i,numRange,dCor3N,dCor3A,alpha,rep);
    power4=calculatePower(power4,i,numRange,dCor4N,dCor4A,alpha,rep);
    power5=calculatePower(power5,i,numRange,dCor5N,dCor5A,alpha,rep);
    power6=calculatePower(power6,i,numRange,dCor6N,dCor6A,alpha,rep);
    power7=calculatePower(power7,i,numRange,dCor7N,dCor7A,alpha,rep);
end

% % Save the results
% filename=strcat('CorrIndTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
% save(filename,'power1','power2','power3','power4','power5','power6','power7','type','n','numRange','rep','lim','dim','noise','option');

function   power1=calculatePower(power1,i,numRange,dCor1N,dCor1A,alpha,rep)
% An auxiliary function to estimate the power based on the distribution of
% the test statistic under the null and the alternative.
k=size(dCor1N,1);
if k>1
    for k1=1:numRange(i);
        for k2=1:numRange(i);
            dCorT=sort(dCor1N(k1,k2,:),'descend');
            cut1=dCorT(ceil(rep*alpha));
            power1(k1,k2,i)=mean(dCor1A(k1,k2,:)>cut1);
        end
    end
else
    dCorT=sort(dCor1N,'descend');
    cut1=dCorT(ceil(rep*alpha));
    power1(i)=mean(dCor1A>cut1);
end
function [power1, power2, power3, power4]=CorrIndTest(type,n,dim,lim,rep, noise)
% Author: Cencheng Shen
% Independence Tests for identifying dependency, returning empirical testing power with respect to increasing sample size at a fixed dimension.
% The output are the powers of local original dCorr, local modified
% dCorr, HHG, and Mantel test.
%
% Parameters:
% type specifies the type of distribution,
% n is the sample size, dim is the dimension
% lim specifies the number of intervals in the sample size,
% rep specifies the number of MC-replicates;
% noise specifies the noise level, by default 0.
if nargin<6
    noise=1; % Default noise level
end
K=n;
alpha=0.05; % Type 1 error level
option1=0; option2=1; option3=1; option4=1; % Control whether to calculate the respective correlation statistic or not.

if lim==0
    numRange=n; % Test at sample size n only
else
    numRange=ceil(n/lim):ceil(n/lim):n; % Test using sample sizes at the interval of ceil(n/lim).
end
lim=length(numRange);

% Intermediate results
dCor1N=zeros(K,K,rep);dCor2N=zeros(K,K,rep);dCor3N=zeros(1,rep);dCor4N=zeros(1,rep);
dCor1A=zeros(K,K,rep);dCor2A=zeros(K,K,rep);dCor3A=zeros(1,rep);dCor4A=zeros(1,rep);
DataN=zeros(n,2*n,rep); DataA=zeros(n,2*n,rep);
d=dim;

% Output
power1=zeros(K,K,lim);power2=zeros(K,K,lim);power3=zeros(1,lim);power4=zeros(1,lim);% Powers for dCorr, mdCorr, HHG, and Mantel

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

for i=1:lim
    nn=numRange(i);
    % First  estimate the distribution of the test statistics under the null
    for r=1:rep
        C=DataN(1:nn,1:nn,r);
        P=DataN(1:nn,n+1:n+nn,r);
        disRank=[disToRanks(C) disToRanks(P)];
        if option1~=0
            dCor1N(1:nn,1:nn,r)=localDCorr(C,P,0,disRank);
        end
        if option2~=0
            dCor2N(1:nn,1:nn,r)=localDCorr(C,P,1,disRank);
        end
        if option3~=0
            dCor3N(r)=HHG(C,P);
        end
        if option4~=0
            dCor4N(r)=Mantel(C,P);
        end
    end
    
    % Then estimate the distribution of the test statistics under the alternative
    for r=1:rep
        C=DataA(1:nn,1:nn,r);
        P=DataA(1:nn,n+1:n+nn,r);
        disRank=[disToRanks(C) disToRanks(P)];
        if option1~=0
            dCor1A(1:nn,1:nn,r)=localDCorr(C,P,0,disRank);
        end
        if option2~=0
            dCor2A(1:nn,1:nn,r)=localDCorr(C,P,1,disRank);
        end
        if option3~=0
            dCor3A(r)=HHG(C,P);
        end
        if option4~=0
            dCor4A(r)=Mantel(C,P);
        end
    end
    
    % Based on the emprical test statistics under the null and the alternative,
    % calculate the emprical powers at type 1 error level alpha=0.05
    for k=1:numRange(i);
        for k2=1:numRange(i);
            dCorT=sort(dCor1N(k,k2,:),'descend');
            cut1=dCorT(ceil(rep*alpha));
            power1(k,k2,i)=mean(dCor1A(k,k2,:)>cut1);
            
            dCorT=sort(dCor2N(k,k2,:),'descend');
            cut2=dCorT(ceil(rep*alpha));
            power2(k,k2,i)=mean(dCor2A(k,k2,:)>cut2);
        end
    end
    dCorT=sort(dCor3N,'descend');
    cut3=dCorT(ceil(rep*alpha));
    power3(i)=mean(dCor3A>cut3);
    dCorT=sort(dCor4N,'descend');
    cut4=dCorT(ceil(rep*alpha));
    power4(i)=mean(dCor4A>cut4);
end

% Save the results
filename=strcat('CorrIndTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
save(filename,'power1','power2','power3','power4','type','n','numRange','rep','lim','dim','noise');

%%% Display and save picture. By default do not display.
% for i=1:length(numRange)
%     power1L(i)=power1(numRange(i),numRange(i),i);
%     power2L(i)=power2(numRange(i),numRange(i),i);
%     power1M(i)=max(max(power1(:,:,i)));
%     power2M(i)=max(max(power2(:,:,i)));
% end
% figure
% plot(numRange,power2M,'ro-',numRange, power1M,'bx-',numRange,power2L,'r.:',numRange,power1L,'b.:',numRange,power3,'g.-',numRange,power4,'c.-','LineWidth',2);
% % Figure title/labels
% hl=legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','SouthEast');
% set(hl,'FontSize',12);
% xlabel('Sample Size','FontSize',12);
% ylabel('Empirical Testing Power','FontSize',12);
% xlim([numRange(1) numRange(end)]);
% ylim([0 1]);
% %titleStr = strcat('Independence Test for ', titlechar,' up to n=', num2str(n), ' at m= ',num2str(dim));
% %title(titleStr,'FontSize',15);
% %saveas(gcf,filename,'jpeg');
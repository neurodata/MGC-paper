function [p1, p2, p3, p4,neighbor1,neighbor2,neighbor3]=CorrPermDistTest(type,rep,cv,titlechar, allP, option)
% Author: Cencheng Shen
% Permutation Tests for identifying dependency.
% The output are the p-values of LGC by mcorr/dcorr/Mantel, and HHG;
% followed by the estimated optimal neighborhood for LGC by
% mcorr/dcorr/Mantel.
% Note that the local family include the global test at the last entry.

% Parameters:
% type should be a n*2n matrix, a concatenation of two distance matrices,
% rep specifies the number of random permutations to use,
% cv specifies the number of bootstrap samples to use for neighborhood validation,
% set allP to non-zero will use all permutations instead,
% option specifies whether each test statistic is calculated or not.
if nargin<3
    cv=1000; % Default bootstrap replicates to estimate the optimal neighborhood 
end
if nargin<4
    titlechar=' Real Data';
end
if nargin<5
    allP=0; % If set to other value, will override rep and use all permutations; unfeasible for n large.
end
if nargin<6
    option=[1,1,1,1]; % Default option. Setting any to 0 to disable the calculation of mcorr/dcorr/Mantel/HHG.
end
n=size(type,1);
C=type(:, 1:n);
P=type(:, n+1:2*n);

% If cv is not 0, use resampling to estimate the optimal neighborhood by
% the testing powers
ps1=zeros(n,n);ps2=zeros(n,n);ps3=zeros(n,n);
if cv~=0
    [p1,p2,p3]=IndependenceTest(C,P,cv);
    neighbor1=verifyNeighbors(1-p1,0);
    neighbor2=verifyNeighbors(1-p2,0);
    neighbor3=verifyNeighbors(1-p3,0);
end
% Return p-values from the permutation test
[p1, p2, p3, p4]=PermutationTest(C,P,rep,allP,option);
if cv==0
    neighbor1=verifyNeighbors(p1,0);
    neighbor2=verifyNeighbors(p2,0);
    neighbor3=verifyNeighbors(p3,0);
end

% Save the results
filename=strcat('CorrPermDistTestType',titlechar);
save(filename,'titlechar','p1','p2','p3','p4','neighbor1','neighbor2','neighbor3','type','n','rep','allP','option');

% %% Plot heatmap
% figure
% K=n;
% kmin=1;
% if n>50
%         c=2;
%         K=ceil(K/2);
%     else
%         c=1;
%         kmin=2;
% end
% xaxis=kmin:K;
% yaxis=kmin:K;
% [X,Y]=meshgrid(c*xaxis,c*yaxis);
% ph=p1(c*xaxis,c*yaxis)';
% surf(X,Y,ph);
% view(2)
% colormap(flipud(colormap))
% caxis([min(min(ph)) 0.1])
% colorbar
% xlabel('Neighborhood Choice of X','FontSize',15);
% ylabel('Neighborhood Choice of Y','FontSize',15);
% xlim([c*kmin,c*K]);
% ylim([c*kmin,c*K]);
% 
% % Figure title/labels
% titleStr = strcat('P-value of MGC for ', titlechar);
% title(titleStr,'FontSize',13);
% filename=strcat('CorrPermDistTest',titlechar);
% saveas(gcf,filename,'jpeg');

function  [p1, p2, p3, p4]=PermutationTest(C,P,rep,allP,option)
% The permutation Test for given data, an auxiliary function of the main
% CorrPermDIstTest function
n=size(C,1);
p1=zeros(n,n); p2=zeros(n,n);p3=zeros(n,n);p4=0;% P-values for LGC by mcorr, LGC by dcorr, LGC by Mantel, and HHG
if nargin<5
    option=[1,1,1,1]; 
end
if allP~=0 
    PAll=perms(1:n);
    rep=size(PAll,1);
end

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
for r2=1:rep
    % Use random permutations; if allP is not 0, use all possible permutations
    per=randperm(n);
    if allP~=0
        per=PAll(r2,:);
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
p1=1-p1;
p2=1-p2;
p3=1-p3;
p4=1-p4;

% Treat the p-value of local methods in neighborhood 0 as 1
% p1(1,:)=1;p1(:,1)=1;p2(1,:)=1;p2(:,1)=1;p3(1,:)=1;p3(:,1)=1;

function  [power1, power2, power3]=IndependenceTest(C,P,rep)
% The independence Test for given data, an auxiliary function of the main
% CorrPermDIstTest function
n=size(C,1);
alpha=0.05; % type 1 error level
% Test statiscs under the null and the alternative
dCor1N=zeros(n,n,rep);dCor2N=zeros(n,n,rep);dCor3N=zeros(n,n,rep);
dCor1A=zeros(n,n,rep);dCor2A=zeros(n,n,rep);dCor3A=zeros(n,n,rep);
% Powers for LGC by mcorr/dcorr/Mantel
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

% Set the powers of all local tests at rank 1 to 0 
power1(1,:)=0;power1(:,1)=0;power2(1,:)=0;power2(:,1)=0;power3(1,:)=0;power3(:,1)=0;
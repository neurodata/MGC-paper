function [p1, p2, p3, p4,neighbor1,neighbor2,neighbor3]=CorrPermDistTest(type,rep,cv,titlechar, allP, option)
% Author: Cencheng Shen
% Permutation Tests for identifying dependency, returning p-value of given data
% The output are the p-values of local original dCorr, local modified dCorr, HHG, and Mantel test.

% Parameters:
% type should be a n*2n matrix, a concatenation of two distance matrices,
% rep specifies the number of random permutations to use,
% cv specifies the number of bootstrap samples to use for neighborhood validation,
% set allP to non-zero will use all permutations instead,
% option specifies whether each test statistic is calculated or not.
if nargin<3
    cv=0; % If set to other value, will override rep and use all permutations; unfeasible for n larger than 8.
end
if nargin<4
    titlechar=' Real Data';
end
if nargin<5
    allP=0; % If set to other value, will override rep and use all permutations; unfeasible for n larger than 8.
end
if nargin<6
    option=[1,1,1,1]; % Default option
end
n=size(type,1);
C=type(:, 1:n);
P=type(:, n+1:2*n);

ps1=zeros(n,n);ps2=zeros(n,n);ps3=zeros(n,n);
if cv~=0
    [p1, p2,p3]=IndependenceTest(C,P,cv);
    neighbor1=verifyNeighbors(p1,0);
    neighbor2=verifyNeighbors(p2,0);
    neighbor3=verifyNeighbors(p3,0);
end
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
% titleStr = strcat('P-value of LGC for ', titlechar);
% title(titleStr,'FontSize',13);
% filename=strcat('CorrPermDistTest',titlechar);
% saveas(gcf,filename,'jpeg');

function  [p1, p2, p3, p4]=PermutationTest(C,P,rep,allP,option)
% Output
n=size(C,1);
p1=zeros(n,n); p2=zeros(n,n);p3=zeros(n,n);p4=0;% P-values for local original dCorr, local modified dCorr, HHG, and Mantel test
if nargin<5
    option=[1,1,1,0];
end
if allP~=0
    PAll=perms(1:n);
    rep=size(PAll,1);
end

% Calculate the test statistics for the given data sets
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

function  [p1, p2, p3]=IndependenceTest(C,P,rep)
% Output
% MDS=0;
% if MDS~=0
%     CEuc=SMDS(C,MDS,0)';
%     PEuc=SMDS(P,MDS,0)';
% end
n=size(C,1);
ratio=1/sqrt(n);
alpha=0.05;
dCor1N=zeros(n,n,rep);dCor2N=zeros(n,n,rep);dCor3N=zeros(n,n,rep);
dCor1A=zeros(n,n,rep);dCor2A=zeros(n,n,rep);dCor3A=zeros(n,n,rep);
power1=zeros(n,n);power2=zeros(n,n);power3=zeros(n,n);
disRankC=disToRanks(C);
for r=1:rep
    per=randsample(n,n,true);
    %per=1:n;
%     per=randperm(n);
%     per2=randperm(n);
%     if MDS==0
%         noise=random('norm',0,1,n,1);
%         noise=squareform(pdist(noise));
%         noise=noise/norm(noise,'fro')*norm(P,'fro')*ratio;
        Pa=P(per,per);%+noise;
        %Pa=(P(per,per)+P(per2,per2))/2;
%         noise=random('norm',0,1,n,1);
%         noise=squareform(pdist(noise));
%         noise=noise/norm(noise,'fro')*norm(C,'fro')*ratio;
        Ca=C(per,per);%+noise;
        %Ca=(C(per,per)+C(per2,per2))/2;
%     else
%         noise=mvnrnd(zeros(n,MDS),eye(MDS,MDS));
%         noise=noise/norm(noise,'fro')*norm(PEuc,'fro')*ratio;
%         Pa=PEuc(per,:);%+noise;
%         %Pa=PEuc(per,:)+PEuc(per2,:);
%         noise=mvnrnd(zeros(n,MDS),eye(MDS,MDS));
%         noise=noise/norm(noise,'fro')*norm(CEuc,'fro')*ratio;
%         Ca=CEuc(per,:);%+noise;
%         %Ca=CEuc(per,:)+CEuc(per2,:);
%         Ca=squareform(pdist(Ca));
%         Pa=squareform(pdist(Pa));
%     end
    disRankC=disToRanks(Ca);
    disRankP=disToRanks(Pa);
    %disRank=[disRankC(per,per) disRankP(per,per)];
    disRank=[disRankC disRankP];
    dCor1A(:,:,r)=LocalGraphCorr(Ca,Pa,1,disRank);
    dCor2A(:,:,r)=LocalGraphCorr(Ca,Pa,2,disRank);
    dCor3A(:,:,r)=LocalGraphCorr(Ca,Pa,3,disRank);
    
    %perN=randsample(n,n,true);
    perN=randperm(n);
%     noise=random('norm',0,1,n,1);
%     noise=squareform(pdist(noise));
%     noise=noise/norm(noise,'fro')*norm(P,'fro')*ratio;
    Pa=P(perN,perN);%+noise;
    disRankP=disToRanks(Pa);
    disRank=[disRankC disRankP];
    dCor1N(:,:,r)=LocalGraphCorr(Ca,Pa,1,disRank);
    dCor2N(:,:,r)=LocalGraphCorr(Ca,Pa,2,disRank);
    dCor3N(:,:,r)=LocalGraphCorr(Ca,Pa,3,disRank);
end

for k=1:n
    for k2=1:n;
        dCorT=sort(dCor1N(k,k2,:),'descend');
        cut1=dCorT(ceil(rep*alpha));
        power1(k,k2)=mean(dCor1A(k,k2,:)>cut1);
        
        dCorT=sort(dCor2N(k,k2,:),'descend');
        cut2=dCorT(ceil(rep*alpha));
        power2(k,k2)=mean(dCor2A(k,k2,:)>cut2);
        
        dCorT=sort(dCor3N(k,k2,:),'descend');
        cut3=dCorT(ceil(rep*alpha));
        power3(k,k2)=mean(dCor3A(k,k2,:)>cut3);
    end
end
power1(1,:)=0;power1(:,1)=0;power2(1,:)=0;power2(:,1)=0;power3(1,:)=0;power3(:,1)=0;
p1=1-power1;
p2=1-power2;
p3=1-power3;

% filename=strcat('CorrPermDistTestTypeN');
% save(filename,'dCor1A','dCor1N','dCor2A','dCor2N','power1','power2');
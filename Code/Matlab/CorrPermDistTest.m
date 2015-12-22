function [p1, p2, p3, p4,neighbor1,neighbor2]=CorrPermDistTest(type,rep,cv,titlechar, allP, option)
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

ps1=zeros(n,n);ps2=zeros(n,n);
if cv~=0
%     %ratio=1/sqrt(n);
%     ratio=1;
%     alpha=0.05;
%     optionCV=option;optionCV(3:4)=0;
%     for i=1:cv
%         noise=squareform(pdist(mvnrnd(zeros(n,n),eye(n,n))));
%         noise=noise/norm(noise,'fro')*norm(P,'fro')*ratio;
%         Pa=P+noise;
%         %per=randsample(n,n,true);
%         per=1:n;
%         %per=randperm(n);
%         [p1, p2,~,~]=PermutationTest(C(per,per),Pa(per,per),rep,allP,optionCV);
%         %ps1=ps1+(p1<alpha)/rep;
%         %ps2=ps2+(p2<alpha)/rep;
%         ps1=ps1+p1/cv;
%         ps2=ps2+p2/cv;
%     end
%     %ps1=1-ps1;
%     %ps2=1-ps2;
    [ps1, ps2]=IndependenceTest(C,P,cv);
end
[p1, p2, p3, p4]=PermutationTest(C,P,rep,allP,option);
if cv==0
    ps1=ps1+p1;
    ps2=ps2+p2;
end
neighbor1=verifyNeighbors(ps1,0);
neighbor2=verifyNeighbors(ps2,0);

% Save the results
filename=strcat('CorrPermDistTestType',titlechar);
save(filename,'titlechar','p1','p2','p3','p4','neighbor1','neighbor2','type','n','rep','allP','option');

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
% surf(X,Y,p2(c*xaxis,c*yaxis));
% view(2)
% colormap(flipud(colormap))
% caxis([0 0.3])
% colorbar
% xlabel('Neighborhood Choice of K','FontSize',15);
% ylabel('Neighborhood Choice of L','FontSize',15);
% xlim([c*kmin,c*K]);
% ylim([c*kmin,c*K]);
% 
% % Figure title/labels
% titleStr = strcat('P-value of Local Graph Dependency for ', titlechar);
% title(titleStr,'FontSize',13);
% filename=strcat('CorrPermDistTest',titlechar);
% saveas(gcf,filename,'jpeg');



% 
% 
% Display and save picture. By default do not display.
% Plot the power/p-value w.r.t. neighborhood
% figure
% K=n;
% kmin=2;xaxis=kmin:K;logOpt=1;
% if logOpt==0
%     plot(xaxis,min(p2(kmin:end,kmin:end),[],2),'r.-',xaxis, min(p1(kmin:end,kmin:end),[],2),'b.-',xaxis,p2(end,end)*ones(length(xaxis),1),'r.:',xaxis,p1(end,end)*ones(length(xaxis),1),'b.:',xaxis,p3*ones(length(xaxis),1),'c.:',xaxis,p4*ones(length(xaxis),1),'g.:','LineWidth',2);
%     ylabel('P-Value','FontSize',15);
%     ylim([0 0.3]);
%     legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','NorthEast');
% else
%     semilogy(xaxis,min(p2(kmin:end,kmin:end),[],2),'r.-',xaxis, min(p1(kmin:end,kmin:end),[],2),'b.-',xaxis,p2(end,end)*ones(length(xaxis),1),'r.:',xaxis,p1(end,end)*ones(length(xaxis),1),'b.:',xaxis,p3*ones(length(xaxis),1),'c.:',xaxis,p4*ones(length(xaxis),1),'g.:','LineWidth',2);
%     ylabel('P-Value','FontSize',15);
%     %ylim([0 1]);
%     legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','NorthEast');
% end
% xlabel('Neighborhood Size','FontSize',15);
% xlim([kmin K]);
% 
% % Figure title/labels
% titleStr = strcat('Permutation Test for ', titlechar);
% title(titleStr,'FontSize',15);
% filename=strcat('CorrPermDistTest',titlechar);
% saveas(gcf,filename,'jpeg');
% alpha=0.05;
% p2(p2>alpha)=1;
% p2=1-p2;
% surf(p2)
% xlim([2,n]);
% ylim([2,n]);
% %title('P-value of Local Modified Distance Correlation','FontSize',15)
% xlabel('Neighborhood k','FontSize',15);
% ylabel('Neighborhood l','FontSize',15);

function  [p1, p2, p3, p4]=PermutationTest(C,P,rep,allP,option)
% Output
n=size(C,1);
p1=zeros(n,n); p2=zeros(n,n);p3=0;p4=0;% P-values for local original dCorr, local modified dCorr, HHG, and Mantel test
if nargin<5
    option=[1,1,0,0];
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
    cut1=localDCorr(C,P,0,disRank);
end
if option(2)~=0
    cut2=localDCorr(C,P,1,disRank);
end
if option(3)~=0
    cut3=HHG(C,P);
end
if option(4)~=0
    cut4=Mantel(C,P);
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
        dCor1=localDCorr(C,Pa,0,disRank);
        p1=p1+(dCor1<cut1)/rep;
    end
    if option(2)~=0
        dCor2=localDCorr(C,Pa,1,disRank);
        p2=p2+(dCor2<cut2)/rep;
    end
    if option(3)~=0
        dCor3=HHG(C,Pa);
        p3=p3+(dCor3<cut3)/rep;
    end
    if option(4)~=0
        dCor4=Mantel(C,Pa);
        p4=p4+(dCor4<cut4)/rep;
    end
end

% Output the p-value
p1=1-p1;
p2=1-p2;
p3=1-p3;
p4=1-p4;

% Treat the p-value of local methods in neighborhood 0 as 1
p1(1,:)=1;p1(:,1)=1;p2(1,:)=1;p2(:,1)=1;

function  [p1, p2]=IndependenceTest(C,P,rep)
% Output
% MDS=0;
% if MDS~=0
%     CEuc=SMDS(C,MDS,0)';
%     PEuc=SMDS(P,MDS,0)';
% end
n=size(C,1);
ratio=1/n;
alpha=0.05;
dCor1N=zeros(n,n,rep);dCor2N=zeros(n,n,rep);
dCor1A=zeros(n,n,rep);dCor2A=zeros(n,n,rep);
power1=zeros(n,n);power2=zeros(n,n);
for r=1:rep
    per=randsample(n,n,true);
    %perN=randsample(n,n,true);
    %per=1:n;
%     per=randperm(n);
%     per2=randperm(n);
    perN=randperm(n);
%     if MDS==0
        noise=random('norm',0,1,n,1);
        noise=squareform(pdist(noise));
        noise=noise/norm(noise,'fro')*norm(P,'fro')*ratio;
        Pa=P(per,per)+noise;
        %Pa=(P(per,per)+P(per2,per2))/2;
        noise=random('norm',0,1,n,1);
        noise=squareform(pdist(noise));
        noise=noise/norm(noise,'fro')*norm(C,'fro')*ratio;
        Ca=C(per,per)+noise;
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
    dCor1A(:,:,r)=localDCorr(Ca,Pa,0,disRank);
    dCor2A(:,:,r)=localDCorr(Ca,Pa,1,disRank);
    disRank=[disRankC disRankP(perN,perN)];
    dCor1N(:,:,r)=localDCorr(Ca,Pa(perN,perN),0,disRank);
    dCor2N(:,:,r)=localDCorr(Ca,Pa(perN,perN),1,disRank);
end

for k=1:n
    for k2=1:n;
        dCorT=sort(dCor1N(k,k2,:),'descend');
        cut1=dCorT(ceil(rep*alpha));
        power1(k,k2)=mean(dCor1A(k,k2,:)>cut1);
        
        dCorT=sort(dCor2N(k,k2,:),'descend');
        cut2=dCorT(ceil(rep*alpha));
        power2(k,k2)=mean(dCor2A(k,k2,:)>cut2);
    end
end
power1(1,:)=0;power1(:,1)=0;power2(1,:)=0;power2(:,1)=0;
p1=1-power1;
p2=1-power2;

% filename=strcat('CorrPermDistTestTypeN');
% save(filename,'dCor1A','dCor1N','dCor2A','dCor2N','power1','power2');
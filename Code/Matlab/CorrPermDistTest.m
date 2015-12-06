function [p1, p2, p3, p4,neighbor1,neighbor2]=CorrPermDistTest(type,rep,titlechar, allP)
% Author: Cencheng Shen
% Permutation Tests for identifying dependency, returning p-value of given data
% The output are the p-values of local original dCorr, local modified
% dCorr, HHG, and Mantel test.

% Parameters:
% type should be a n*2n matrix, a concatenation of two distance matrices,
% rep specifies the number of random permutations to use,
% set allP to non-zero will use all permutations instead.
if nargin<3
    titlechar=' Real Data';
end
if nargin<4
    allP=0; % If set to other value, will override rep and use all permutations; unfeasible for n larger than 8.
end
n=size(type,1);
K=n;
%alpha=0.05; %type 1 error level
option1=1; option2=1; option3=1; option4=1; % Control whether to calculate the respective correlation statistic or not.
cut1=0;cut2=0;cut3=0;cut4=0;

% Output
p1=zeros(K,K); p2=zeros(K,K);p3=0;p4=0;% P-values for local original dCorr, local modified dCorr, HHG, and Mantel test

% Calculate the test statistics for the given data sets
C=type(:, 1:n);
P=type(:, n+1:2*n);
disRankC=disToRanks(C);
disRankP=disToRanks(P);
disRank=[disRankC disRankP];
if option1~=0
    cut1=localDCorr(C,P,0,disRank);
end
if option2~=0
    cut2=localDCorr(C,P,1,disRank);
end
if option3~=0
    cut3=HHG(C,P);
end
if option4~=0
    cut4=Mantel(C,P);
end

% Permute the second dataset for rep times, and calculate all test
% statistics between the permuted data sets
if allP~=0
    PAll=perms(1:n);
    rep=size(PAll,1);
end
for r2=1:rep
    % Use random permutations; if allP is not 0, use all possible permutations
    per=randperm(n);
    if allP~=0
        per=PAll(r2,:);
    end
    Pa=P(per,per);
    disRank=[disRankC disRankP(per, per)];
    if option1~=0
        dCor1=localDCorr(C,Pa,0,disRank);
        p1=p1+(dCor1<cut1)/rep;
    end
    if option2~=0
        dCor2=localDCorr(C,Pa,1,disRank);
        p2=p2+(dCor2<cut2)/rep;
    end
    if option3~=0
        dCor3=HHG(C,Pa);
        p3=p3+(dCor3<cut3)/rep;
    end
    if option4~=0
        dCor4=Mantel(C,Pa);
        p4=p4+(dCor4<cut4)/rep;
    end
end
% Find the best neighborhood
%if cv==0
    neighbor1=findNeighbor(p1);
    neighbor2=findNeighbor(p2);
%end

% Output the p-value
p1=1-p1;
p2=1-p2;
p3=1-p3;
p4=1-p4;

% Save the results
filename=strcat('CorrPermDistTestType',titlechar);
save(filename,'titlechar','p1','p2','p3','p4','neighbor1','neighbor2','cut1','cut2','cut3','cut4','type','n','rep');

% Display and save picture. By default do not display.
%%Plot the power/p-value w.r.t. neighborhood
% p1(neighbor1(1), neighbor1(2))
% p2(neighbor2(1), neighbor2(2))
% p1(end,end)
% p2(end,end)
% p3
% p4
% figure
% K=n;
% kmin=1;xaxis=kmin:K;logOpt=1;
% if logOpt==0
%     plot(xaxis,min(p2,[],2),'r.-',xaxis, min(p1,[],2),'b.-',xaxis,p2(end,end)*ones(length(xaxis),1),'r.:',xaxis,p1(end,end)*ones(length(xaxis),1),'b.:',xaxis,p3*ones(length(xaxis),1),'c.-',xaxis,p4*ones(length(xaxis),1),'g.-','LineWidth',2);
%     ylabel('P-Value','FontSize',15);
%     ylim([0 1]);
%     legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','NorthEast');
% else
%     plot(xaxis,log(min(p2,[],2)),'ro-',xaxis, log(min(p1,[],2)),'bx-',xaxis,log(p2(end,end))*ones(length(xaxis),1),'r.:',xaxis, log(p1(end,end))*ones(length(xaxis),1),'b.:',xaxis,log(p3)*ones(length(xaxis),1),'g.-',xaxis,log(p4)*ones(length(xaxis),1),'c.-','LineWidth',2);
%     ylabel('Log p-value','FontSize',15);
%     ylim([-5 0]);
%     legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','SouthEast');
% end
% xlabel('Neighborhood Size','FontSize',15);
% xlim([kmin K]);
% 
% % Figure title/labels
% titleStr = strcat('Permutation Test for ', titlechar,' at n=', num2str(n));
% title(titleStr,'FontSize',15);
% filename=strcat('CorrPermDistTest',titlechar);
%saveas(gcf,filename,'jpeg');

function neighbor=findNeighbor(p)
neighbor=zeros(2,1);
pmax=max(max(p));
neighbor(2)=find(max(p)==pmax,1);
neighbor(1)=find(p(:,neighbor(2))==pmax,1);

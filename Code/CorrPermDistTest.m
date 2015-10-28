function [p1, p2, p3, p4]=CorrPermDistTest(type,n,lim,rep1,rep2,titlechar)
% Author: Cencheng Shen
% Permutation Tests for identifying correlation, returning p-value with respect to increasing sample size.
% The output is the p-value of rankDCorr, dCorr, modified dCorr, HHG.

% Parameters:
% In the input, rep1 specifies the number of MC replicates.
% rep2 specifies the number of random permutations.
if nargin<6
    titlechar=' Real Data';
end
display=0;
K=n;
alpha=0.05; %type 1 error level
option1=1; option2=1; option3=1; option4=1; %Control whether to calculate the respective correlation statistic or not.

if lim==0
    numRange=n; %do sample size at n
else
    numRange=ceil(n/lim):ceil(n/lim):n; %do sample size at the interval of ceil(n/lim).
end
lim=length(numRange);

% Output
p1=zeros(lim,K,rep1); p2=zeros(lim,K,rep1);p3=zeros(lim,K,rep1);p4=zeros(lim,rep1);%p-values for rankdCorr, dCorr, modified dCorr, HHG
power1=zeros(lim,K,1); power2=zeros(lim,K,1);power3=zeros(lim,K,1);power4=zeros(lim,1);%average powers for rankdCorr, dCorr, modified dCorr, HHG
dCor1=zeros(rep2,K);dCor2=zeros(rep2,K);dCor3=zeros(rep2,K);dCor4=zeros(rep2,1);

%rep1 MC replicates
for r1=1:rep1    
    % Iterate through all sample size choices
    C1=type(:, 1:n);
    P1=type(:, n+1:2*n);
%     C1=(C1+C1')/2;
%     P1=(P1+P1')/2;
    for i=1:lim
        nn=numRange(i);
        C=C1(1:nn,1:nn);
        P=P1(1:nn,1:nn);
        disRankC=disToRanks(C);
        disRankP=disToRanks(P);
        
        % Permute the second dataset for rep2 times, and calculate the 4 correlation statistics
        for r2=1:rep2
            per=randperm(nn);
            Pa=P(per,per);
            disRank=[disRankC disRankP(per, per)];
            if option1~=0
                dCor1(r2,1:nn)=rankDCorr(C,Pa,disRank,1);
            end
            if option2~=0
                dCor2(r2,1:nn)=rankDCorr(C,Pa,disRank,2);
            end
            if option3~=0
                dCor3(r2,1:nn)=rankDCorr(C,Pa,disRank,3);
            end
            if option4~=0
                dCor4(r2)=HHG(C,Pa);
            end
        end
        
        % Calculate the correlation statistics for the original data, and
        % derive the p-value for the permutation test.
        disRank=[disRankC disRankP];
        if option1~=0
            cut1=rankDCorr(C,P,disRank,1);
            p1(i,1:nn,r1)=mean(dCor1(:,1:nn)<repmat(cut1,1,rep2)',1);
        end
        if option2~=0
            cut2=rankDCorr(C,P,disRank,2);
            p2(i,1:nn,r1)=mean(dCor2(:,1:nn)<repmat(cut2,1,rep2)',1);
        end
        if option3~=0
            cut3=rankDCorr(C,P,disRank,3);
            p3(i,1:nn,r1)=mean(dCor3(:,1:nn)<repmat(cut3,1,rep2)',1);
        end
        if option4~=0
            cut4=HHG(C,P);
            p4(i,r1)=mean(dCor4<cut4);
        end
    end
    power1=power1+(p1(:,:,r1)>(1-alpha))/rep1;
    power2=power2+(p2(:,:,r1)>(1-alpha))/rep1;
    power3=power3+(p3(:,:,r1)>(1-alpha))/rep1;
    power4=power4+(p4(:,r1)>(1-alpha))/rep1;
end

% Output the mean p-value and standard deviation. Std is not meaningful
% when rep1=1.
std1=std(p1,0,3);
std2=std(p2,0,3);
std3=std(p3,0,3);
std4=std(p4,0,2);
p1=1-mean(p1,3);
p2=1-mean(p2,3);
p3=1-mean(p3,3);
p4=1-mean(p4,2);

% Save the results
filename=strcat('CorrPermDistTestType',titlechar);
save(filename,'titlechar','power1','power2','power3','power4','p1','p2','p3','p4','std1','std2','std3','std4','dCor1','dCor2','dCor3','dCor4','cut1','cut2','cut3','cut4','type','n','numRange','lim','rep1','rep2');
%save(filename,'titlechar','power1','power2','power3','power4','p1','p2','p3','p4','std1','std2','std3','std4','type','n','numRange','lim','rep1','rep2');

% Display and save picture. By default do not display.
if display~=0
%     figure
%     %%% Plot the average p-value
%     % plot(numRange,min(p1(:,:),[],2),'ro-',numRange, min(p2(:,:),[],2),'bx-',numRange,min(p3(:,:),[],2),'g^-',numRange,p4,'c.-','LineWidth',2);
%     % plot(numRange,min(p1(:,:),[],2),'ro-',numRange, p2(:,end),'bx-',numRange, p3(:,end),'g^-',numRange,p4,'c.-','LineWidth',2);
%     %legend('Rank Distance Correlation','Original Distance Correlation','Modified Distance Correlation','HHG','Location','NorthEast');
%     % Plot the permutation test power
%     plot(numRange,max(power1,[],2),'ro-',numRange, max(power2,[],2),'bx-',numRange,max(power3,[],2),'g^-',numRange,power4,'c.-','LineWidth',2);
%     for i=1:length(numRange)
%         power2L(i)=power2(i,numRange(i)-1);
%         power3L(i)=power3(i,numRange(i)-1);
%     end
%     plot(numRange,max(power1,[],2),'ro-',numRange, power2L,'bx-',numRange, power3L,'g^-',numRange,power4,'c.-','LineWidth',2);
%     legend('Rank Distance Correlation','Original Distance Correlation','Modified Distance Correlation','HHG','Location','SouthEast');
%     xlim([numRange(1) numRange(end)]);
%     ylim([0 1]);
%     xlabel('Sample Size');
%     ylabel('Power');
%     titleStr = strcat('Permutation Test for ', titlechar,' upto n=', num2str(n));
%     title(titleStr);
%     saveas(gcf,filename,'jpeg');
    
    %%Plot the power/p-value w.r.t. neighborhood
    figure
    K=n;
    kmin=1;xaxis=kmin:K;
    plot(xaxis,p1(lim,xaxis),'ro-',xaxis, p2(lim,xaxis),'bx-',xaxis,p3(lim,xaxis),'g^-',xaxis,p4(lim)*ones(length(xaxis),1),'c.-','LineWidth',1);
    %plot(xaxis,p1(lim,xaxis),'ro-',xaxis, p2(lim,end)*ones(length(xaxis),1),'bx-',xaxis,p3(lim,end)*ones(length(xaxis),1),'g^-',xaxis,p4(lim)*ones(length(xaxis),1),'c.-','LineWidth',1);
    xlabel('Neighborhood Size');
    xlim([kmin K]);
    ylabel('P-Value');
    ylim([0 1]);
    
    % Figure title/labels
    legend('Rank Distance Correlation*','Original Distance Correlation*','Modified Distance Correlation*','HHG','Location','NorthEast');
    titleStr = strcat('Permutation Test for ', titlechar,' at n=', num2str(n));
    title(titleStr);
    filename=strcat('CorrPermDistTest',titlechar);
    saveas(gcf,filename,'jpeg');
end
function [p,pmean,pstd]=CORRDataSummary()

fileDir='../../../../../Data/CORR/';
pre1='../Data/'; % The folder to locate data
pre2='../../../Draft/Figures/FigReal'; % The folder to save figures

allFiles = dir( fileDir );
allNames = { allFiles.name };
strPost='FalseDetection';
indStr = strfind(allNames,strPost);
indStr = find(~cellfun(@isempty,indStr));
allNames=allNames(indStr);
n=length(allNames);

strNames=cell(n,1);
p=zeros(n,7);
for j=1:n
    fileName=allNames(j);
    fileName=fileName{1,1};
    strInd=findstr(fileName,strPost);
    strNames{j}=fileName(1:strInd(1)-1);
    load(strcat(fileDir,fileName),'power');
    p(j,:)=power;
end
pmean=mean(p);
pstd=std(p);

%%%sort
[~,ind]=sort(p(:,2),'ascend');

cmap=zeros(3,3);
ma = [1,0,1];
cmap(1,:) = ma;
cmap(2,:) = ma;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

x=1:n;
% hold on
%semilogy(x,p(ind,2),'ob:',x,p(ind,5),'xr:',x,ones(n,1)*pmean(2),'b-',x,ones(n,1)*pmean(5),'r-'); %Jittered data
plot(x,p(ind,2),'.-',x,ones(n,1)*pmean(2),'.:',x,ones(n,1)*0.05,'k.:','LineWidth',2); %Jittered data
xlim([1,n]);
ylim([0,0.3]);
set(gca,'FontSize',14);
ax=gca;
h=legend('MGC','Mean FPR of MGC', 'FPR=0.05', 'Location','NorthEast');
set(h,'FontSize',16);

set(gca,'XTickLabel',strNames(ind),'XTick',1:numel(strNames),'FontSize',11);
ax.XTickLabelRotation=45;
%ylim([0,1]);
ylabel('False Positive Rate','FontSize',16);
title('False Positive Rates for Brain vs Noise','FontSize',17);

F.fname=strcat(pre2, 'CORR');
F.wh=[3 2.5]*2;
print_fig(gcf,F)
save(strcat(pre1,'CORRSummary.mat'),'p','pmean','pstd','n');
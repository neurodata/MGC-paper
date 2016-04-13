function [p,pmean,pstd]=CORRDataSummary()

fileDir='../../../../Data/CORR/';
pre1='../../Data/'; % The folder to locate data
pre2='../../Figures/FigReal'; % The folder to save figures

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

x=1:n;
plot(x,p(:,2),'ob--',x,p(:,5),'xr--'); %Jittered data
xlim([1,n]);
ax=gca;
h=legend('MGC\{mcorr\}','mcorr','Location','NorthEast');
set(h,'FontSize',13);
% strNames=['J';'F';'M';'A';'M'];
set(gca,'XTickLabel',strNames,'XTick',1:numel(strNames),'FontSize',12);
ax.XTickLabelRotation=45;
%xlim([0.8,2.2])
ylim([0,1]);
%xlabel('Data Name','FontSize',13);
ylabel('False Positive Rate','FontSize',13);
% set(gca,'XTick',[]); % Remove x axis ticks
title('MGC False Positive Rate Scatter Plot for Brain vs Noise','FontSize',13);

F.fname=strcat(pre2, num2str(2));
F.wh=[3 2.5]*2;
print_fig(gcf,F)
save(strcat(pre1,'CORRSummary.mat'),'p','pmean','pstd','n');
function [p,pmean,pstd]=CORRDataSummary()

fileDir='../../../../Data/CORR/';
pre1='../../Data/'; % The folder to locate data
pre2='../../Figures/FigReal'; % The folder to save figures

allFiles = dir( fileDir );
allNames = { allFiles.name };
indStr = strfind(allNames,'FalseDetection');
indStr = find(~cellfun(@isempty,indStr));
allNames=allNames(indStr);
n=length(allNames);

p=zeros(n,7);
for j=1:n
    fileName=allNames(j);
    load(strcat(fileDir,fileName{1,1}),'power');
    p(j,:)=power;
end
pmean=mean(p);
pstd=std(p);

x=ones(n,1);
x=x+rand(size(x));
plot(x,p(:,1),'ob') %Jittered data
xlim([0.8,2.2])
ylim([0,0.2])
set(gca,'XTick',[]); % Remove x axis ticks
title('MGC False Positive Rate Scatter Plot for Brain vs Noise','FontSize',13);

F.fname=strcat(pre2, num2str(2));
F.wh=[3 2.5]*2;
print_fig(gcf,F)
save(strcat(pre1,'CORRSummary.mat'),'p','pmean','pstd','n');
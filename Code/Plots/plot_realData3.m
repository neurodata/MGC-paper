function []=plot_realData3
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
pre2=strcat(rootDir,'Figures/FigReal');% The folder to save figures
sameBar=1;

%% figure stuff

cmap=zeros(4,3);
gr =[0,1,0];
ma = [1,0,1];
map3 = brewermap(128,'PiYG'); % brewmap
lgr=map3(100,:);
dgr=map3(128,:);
gr=lgr;
glob= [0.5,0.5,0.5];
% cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = gr;
cmap(3,:) = gr;
cmap(4,:) = gr;
% cmap(3,:) = cy;
map1=cmap;
map2 = brewermap(128,'BuPu'); % brewmap
set(groot,'defaultAxesColorOrder',map1);

%% Mantel

type=1;
    load(strcat(rootDir,'Data/Preprocessed/BrainCP.mat'))
    load(strcat(rootDir,'Data/Results/CorrPermDistTestTypeBrainCxP.mat'))
    n=42;
    C=distC;
    D=distP;
    fontSize=15;
pos1 =[0.05,0.6,0.33,0.33];
pos2 =[0.55,0.6,0.33,0.33];
pos3 =[0.05,0.1,0.33,0.33];
pos4 =[0.55,0.1,0.33,0.33];
ax=subplot('position',pos1);
hold on
kmin=1;
[testMGC, testMLocal, optimalInd]=MGCSampleStat(C,D);
[m,n]=size(testMLocal);
[k,l]=ind2sub(size(testMLocal),optimalInd(1));
ph=testMLocal(kmin:m,kmin:n)';
%indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% ph(indPower)=2;
imagesc(ph);
plot(m,n,'.','markerSize',40,'MarkerFaceColor',glob,'Color',glob)
plot(k,l,'g.','markerSize',40)
hold off
axis('square')
pos = get(ax,'position');

set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
% cmap=map4;
colormap(ax,map2)
hm=ceil(max(max(ph))*100)/100;
% hm=ceil(prctile(ph(ph<1),99)*100)/100;
caxis([0 hm])
h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
set(h,'FontSize',fontSize);
% h2=get(h,'Position');
% h2(4)=h2(4)/2;
% set(h,'Position',h2)
xlim([1 m]);
ylim([1 n]);
%set(gca,'XTick',[2.5,round(n/2)-1,n],'YTick',[2.5,round(n/2)-1,n],'XTickLabel',[2,round(m/2),n],'YTickLabel',[2,round(m/2),n],'FontSize',16);
%set(gca,'XTick',[],'YTick',[])
% pos = get(ax,'position');
title('i. Brain vs Personality','FontSize',fontSize-1, 'Units', 'normalized', ...
    'Position', [0 1.05], 'HorizontalAlignment', 'left');
% title('Local Correlations','fontweight','normal','FontSize',fontSize);
xlabel('# X Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
    'Position', [0 -0.17], 'HorizontalAlignment', 'left')
ylabel('# Y Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
    'Position', [-0.2 0.25], 'VerticalAlignment', 'bottom')
%text(10,110,'Test Statistics','FontSize',fontSize)
% pos2 = get(ax,'position');
% pos2(3:4) = F.pos(3:4);
set(ax,'position',pos);

    load(strcat(rootDir,'Data/Preprocessed/semipar.mat'))
    load(strcat(rootDir,'Data/Results/CorrPermDistTestTypeMigrainxCCI.mat'))
    n=109;
    distCCI=squareform(pdist(cci));
    C=distMigrain(ind,ind);
    D=distCCI(ind,ind);
ax=subplot('position',pos2);
hold on
kmin=1;
[testMGC, testMLocal, optimalInd]=MGCSampleStat(C,D);
[m,n]=size(testMLocal);
[k,l]=ind2sub(size(testMLocal),optimalInd(1));
ph=testMLocal(kmin:m,kmin:n)';
%indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% ph(indPower)=2;
imagesc(ph);
plot(m,n,'.','markerSize',40,'MarkerFaceColor',glob,'Color',glob)
plot(k,l,'g.','markerSize',40)
hold off
axis('square')
pos = get(ax,'position');

set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
% cmap=map4;
colormap(ax,map2)
hm=ceil(max(max(ph))*100)/100;
hm=hm*1;
% hm=ceil(prctile(ph(ph<1),99)*100)/100;
caxis([0 hm])
h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
set(h,'FontSize',fontSize);
% h2=get(h,'Position');
% h2(4)=h2(4)/2;
% set(h,'Position',h2)
xlim([1 m]);
ylim([1 n]);
%set(gca,'XTick',[2.5,round(n/2)-1,n],'YTick',[2.5,round(n/2)-1,n],'XTickLabel',[2,round(m/2),n],'YTickLabel',[2,round(m/2),n],'FontSize',16);
%set(gca,'XTick',[],'YTick',[])
% pos = get(ax,'position');
title('ii. Brain vs Creativity','FontSize',fontSize-1, 'Units', 'normalized', ...
    'Position', [0 1.05], 'HorizontalAlignment', 'left');
% title('Local Correlations','fontweight','normal','FontSize',fontSize);
xlabel('# X Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
    'Position', [0 -0.17], 'HorizontalAlignment', 'left')
ylabel('# Y Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
    'Position', [-0.2 0.25], 'VerticalAlignment', 'bottom')
%text(10,110,'Test Statistics','FontSize',fontSize)
% pos2 = get(ax,'position');
% pos2(3:4) = F.pos(3:4);
set(ax,'position',pos);

%     load(strcat(rootDir,'Data/Preprocessed/proteomics.mat'))
% %     per=(LabelIndAll~=2) & (LabelIndAll<5);
% %     LabelIndAll(per)=1;
%     per=(LabelIndAll<5);
%     m=318;
%     D=LabelIndAll(per);
%     C=A(:,per)';
%     D=squareform(pdist(D));
%     C=squareform(pdist(C));
%     
%     fontSize=15;
% pos1 =[0,0.33,0.33,0.33];
% pos2 =[0.34,0.33,0.33,0.33];
% pos3 =[0.67,0.33,0.33,0.33];
% ax=subplot('position',pos2);
% hold on
% kmin=1;
% [testMGC, testMLocal, optimalInd]=MGCSampleStat(C,D);
% [m,n]=size(testMLocal);
% [k,l]=ind2sub(size(testMGC),optimalInd(1));
% ph=testMLocal(kmin:m,kmin:n)';
% %indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% % ph(indPower)=2;
% imagesc(ph);
% plot(m,n,'.','markerSize',40,'MarkerFaceColor',glob,'Color',glob)
% plot(l,k,'g.','markerSize',40)
% hold off
% axis('square')
% pos = get(ax,'position');
% 
% set(gca,'FontSize',fontSize)
% set(gca,'YDir','normal')
% % cmap=map4;
% colormap(ax,map2)
% hm=ceil(max(max(ph))*100)/100;
% % hm=ceil(prctile(ph(ph<1),99)*100)/100;
% caxis([0 hm])
% h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
% set(h,'FontSize',fontSize);
% % h2=get(h,'Position');
% % h2(4)=h2(4)/2;
% % set(h,'Position',h2)
% xlim([1 m]);
% ylim([1 n]);
% %set(gca,'XTick',[2.5,round(n/2)-1,n],'YTick',[2.5,round(n/2)-1,n],'XTickLabel',[2,round(m/2),n],'YTickLabel',[2,round(m/2),n],'FontSize',16);
% %set(gca,'XTick',[],'YTick',[])
% % pos = get(ax,'position');
% title('MGC Image','FontSize',fontSize-1, 'Units', 'normalized', ...
%     'Position', [0 1.05], 'HorizontalAlignment', 'left');
% % title('Local Correlations','fontweight','normal','FontSize',fontSize);
% xlabel('# X Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
%     'Position', [0 -0.17], 'HorizontalAlignment', 'left')
% ylabel('# Y Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
%     'Position', [-0.15 0.3], 'VerticalAlignment', 'bottom')
% %text(10,110,'Test Statistics','FontSize',fontSize)
% % pos2 = get(ax,'position');
% % pos2(3:4) = F.pos(3:4);
% set(ax,'position',pos);

load(strcat(pre1,'ScreeningPancvsNormal.mat'))
testMGC2=testMGC;
load(strcat(pre1,'ScreeningPancvsAll.mat'))
ax=subplot('position',pos3);
set(ax,'FontSize',fontSize-1);
% [~,ind]=sort(testMGC,'descend');
% plot(testMGC(ind),'.-','LineWidth',lw);
hold on
p1=testMGC2(:,6);
p2=testMGC(:,6);
plot(p1,p2,'.','Color',gr,'MarkerSize',15);
x=0:0.0001:1;
y=0.001;
plot(x,y*ones(length(x),1),'--','Color',glob);
plot(y*ones(length(x),1),x,'--','Color',glob);
ind=181;
plot(p1(ind),p2(ind),'.','Color',gr,'MarkerSize',25);
text(p1(ind),p2(ind),'neurogranin','VerticalAlignment','bottom','HorizontalAlignment','left','Color',gr,'FontSize',fontSize+5);
hold off
set(gca,'XScale','log','YScale','log');
% xlim([1,318]);
% ylim([-0.1,0.5]);
%set(gca,'YTick',[0,0.2,0.4],'XTick',[1,100,200,300]);
% ylabel('Magnitude','FontSize',fs, ...
%     'Units', 'normalized','Position', [-0.21 0], 'HorizontalAlignment', 'left')
% xlabel('Features','FontSize',fs,...
%     'Units', 'normalized','Position', [-0.01, -0.16], 'HorizontalAlignment', 'left')
ylabel('P-values for Panc vs All','FontSize',fontSize, ...
    'Units', 'normalized','Position', [-0.2 0], 'HorizontalAlignment', 'left')
xlabel('P-values for Panc vs Norm','FontSize',fontSize,...
    'Units', 'normalized','Position', [-0.01, -0.17], 'HorizontalAlignment', 'left')
title('iii. Cancer Biomarker Discovery','FontSize',fontSize, ...
    'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left');
axis('square');


ax=subplot('position',pos4);
set(ax,'FontSize',fontSize-1);
% [~,ind]=sort(testMGC,'descend');
% plot(testMGC(ind),'.-','LineWidth',lw);
hold on
y=[2 5; 4 5; 2 1]; %all: 1/3
% cat = categorical({'MGC','HSIC','HHG'});
h=bar(y);
set(h(1),'FaceColor','r') ;
set(h(2),'FaceColor','g') ;
% for k = 1:size(y,2)
%     h(1).CData = 'r';
%     h(2).CData = 'g';
% end
% p1=testMGC2(:,6);
% p2=testMGC(:,6);
% plot(p1,p2,'.','Color',gr,'MarkerSize',15);
% x=0:0.0001:1;
% y=0.001;
% plot(x,y*ones(length(x),1),'--','Color',glob);
% plot(y*ones(length(x),1),x,'--','Color',glob);
% ind=181;
% plot(p1(ind),p2(ind),'.','Color',gr,'MarkerSize',25);
% text(p1(ind),p2(ind),'neurogranin','VerticalAlignment','bottom','HorizontalAlignment','left','Color',gr,'FontSize',fontSize+5);
hold off
% set(gca,'XScale','log','YScale','log');
xlim([0.5,3.5]);
ylim([0,8]);
names = {'MGC'; 'HSIC'; 'HHG'};
set(gca,'XTick',[1:3],'XTickLabel',names);
legend('False Positives','True Positives');
% ylabel('Magnitude','FontSize',fs, ...
%     'Units', 'normalized','Position', [-0.21 0], 'HorizontalAlignment', 'left')
% xlabel('Features','FontSize',fs,...
%     'Units', 'normalized','Position', [-0.01, -0.16], 'HorizontalAlignment', 'left')
ylabel('# True / False Positives','FontSize',fontSize, ...
    'Units', 'normalized','Position', [-0.2 0], 'HorizontalAlignment', 'left')
% xlabel('# of True / False Positives','FontSize',fontSize,...
%      'Units', 'normalized','Position', [-0.01, -0.17], 'HorizontalAlignment', 'left')
title('iv. Biomarker kNN Classification','FontSize',fontSize, ...
    'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left');
axis('square');



pre2=strcat(rootDir,'Figures/');% The folder to save figures
F.fname=strcat(pre2, 'FigReal');
F.wh=[5 4]*2;
F.PaperPositionMode='auto';

print_fig(gcf,F)
%clean_panel(ax,map4,pos,id,n,col,fontSize)
%%
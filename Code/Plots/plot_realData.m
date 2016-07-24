function []=plot_realData
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
pre2=strcat(rootDir,'Figures/FigReal');% The folder to save figures

cmap=zeros(4,3);
gr =[0,1,0];
ma = [1,0,1];
map3 = brewermap(128,'PiYG'); % brewmap
lgr=map3(100,:);
dgr=map3(128,:);
gr=lgr;
% cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = gr;
cmap(3,:) = gr;
cmap(4,:) = gr;
% cmap(3,:) = cy;
map1=cmap;
map2 = brewermap(128,'GnBu'); % brewmap


% Plot 1st figure
set(groot,'defaultAxesColorOrder',map1);
filename=strcat(pre1,'CorrPermDistTestTypeBrainCxP.mat');
load(filename);
n=size(p2All,1);

figure
kmin=2;
imagesc(log(p2All(kmin:end,kmin:end)'));
set(gca,'YDir','normal')
colormap(flipud(map2))
cticks=[0.01, 0.1, 0.5];
h=colorbar('Ticks',log(cticks),'TickLabels',cticks);%,'location','westoutside');
set(gca,'FontSize',24);
%     if i==3
[n1,n2]=size(p2All);
xlabel('# of Neighbors for X','FontSize',24);
ylabel('# of Neighbors for Y','FontSize',24);
set(gca,'XTick',[1,round(n1/2)-1,n1-1],'XTickLabel',[2,round(n1/2),n1]); % Remove x axis ticks
set(gca,'YTick',[1,round(n2/2)-1,n2-1],'YTickLabel',[2,round(n2/2),n2]); % Remove x axis ticks
%h=colorbar('Ticks',[0,0.05,0.1]);%,'location','westoutside');
%     else
%         set(gca,'XTick',[]); % Remove x axis ticks
%         set(gca,'YTick',[]); % Remove y axis ticks
%     end
title('Brain Activity vs. Personality','FontSize',24);

F.fname=strcat(pre2, num2str(1));
F.wh=[3 2.5]*2;
print_fig(gcf,F)

% kmin=2;
% xa=kmin:n;
pp1=p2All(:,end);
filename=strcat(pre1,'CorrPermDistTestTypeBrainLMRxY.mat');
load(filename);
n=size(p2All,1);
pp2=p2All(:,end);

figure
kmin=2;
imagesc(log(p2All(kmin:end,kmin:end)'));
set(gca,'YDir','normal')
colormap(flipud(map2))
cticks=[0.01, 0.1, 0.5];
h=colorbar('Ticks',log(cticks),'TickLabels',cticks);%,'location','westoutside');
set(gca,'FontSize',24);
%     if i==3
[n1,n2]=size(p2All);
xlabel('# of Neighbors for X','FontSize',24);
ylabel('# of Neighbors for Y','FontSize',24);
set(gca,'XTick',[1,round(n1/2)-1,n1-1],'XTickLabel',[2,round(n1/2),n1]); % Remove x axis ticks
set(gca,'YTick',[1,round(n2/2)-1,n2-1],'YTickLabel',[2,round(n2/2),n2]); % Remove x axis ticks
%h=colorbar('Ticks',[0,0.05,0.1]);%,'location','westoutside');
%     else
%         set(gca,'XTick',[]); % Remove x axis ticks
%         set(gca,'YTick',[]); % Remove y axis ticks
%     end
title('Brain Shape vs. Disorder','FontSize',24);

F.fname=strcat(pre2, num2str(2));
F.wh=[3 2.5]*2;
print_fig(gcf,F)

% figure
% y=0.05;
% H=area(xa,y*ones(n-kmin+1,1),'LineStyle',':','FaceColor', [.9 .9 .9]);
% hold on
% plot(xa,pp1,'k.-',xa,pp2,'.--','LineWidth',2)
% set(gca,'FontSize',18);
% xlabel('Number of Neighbors for X');
% xlim([2 n]);
% ylim([0 0.2]);
% ylabel('P-Value');

% figure
% hold on
% ind1=1:n;
% ind1=ind1(pp1>0.05);
% ind1=[1 ind1 n+1];
% [~,i]=max(ind1(2:end)-ind1(1:length(ind1)-1));
% ind1=ind1(i)+1:ind1(i+1)-2;
% ind1S=2:ind1(1);
% ind1E=ind1(end):n;
% 
% ind2=1:n;
% ind2=ind2(pp2>0.05);
% ind2=[1 ind2 n+1];
% [~,i]=max(ind2(2:end)-ind2(1:length(ind2)-1));
% ind2=ind2(i)+1:ind2(i+1)-2;
% ind2S=2:ind2(1);
% ind2E=ind2(end):n;
% 
% plot(ind1S,pp1(ind1S),'-',ind2S,pp2(ind2S),':',ind1E,pp1(ind1E),'-',ind2E,pp2(ind2E),':','LineWidth',3,'Color',lgr);
% h=legend('Left Brain','Right Brain','Location','North');
% set(h,'FontSize',24);
% legend boxoff
% plot(ind1,pp1(ind1),'k-','Color',dgr,'LineWidth',6)
% plot(ind2,pp2(ind2),'k:','Color',dgr,'LineWidth',6)
% set(gca,'FontSize',24);
% xlabel('# of Neighbors for X','FontSize',24);
% n1=length(pp1);
% set(gca,'XTick',[2,round(n1/2),n1],'XTickLabel',[2,round(n1/2),n1]); % Remove x axis ticks
% set(gca,'YTick',[0 0.05 0.1 0.15 0.2]); % Remove x axis ticks
% xlim([2 n]);
% ylim([0 0.2]);
% ylabel('P-Value','FontSize',24);

% Plot second figure
filename=strcat(pre1,'CorrPermDistTestTypeMigrainxCCI.mat');
load(filename);
figure
kmin=2;
imagesc(log(p2All(kmin:end,kmin:end)'));
set(gca,'YDir','normal')
colormap(flipud(map2))
cticks=[0.01, 0.1, 0.5];
h=colorbar('Ticks',log(cticks),'TickLabels',cticks);%,'location','westoutside');
set(gca,'FontSize',24);
%     if i==3
[n1,n2]=size(p2All);
xlabel('# of Neighbors for X','FontSize',24);
ylabel('# of Neighbors for Y','FontSize',24);
set(gca,'XTick',[1,round(n1/2)-1,n1-1],'XTickLabel',[2,round(n1/2),n1]); % Remove x axis ticks
set(gca,'YTick',[1,round(n2/2)-1,n2-1],'YTickLabel',[2,round(n2/2),n2]); % Remove x axis ticks
%h=colorbar('Ticks',[0,0.05,0.1]);%,'location','westoutside');
%     else
%         set(gca,'XTick',[]); % Remove x axis ticks
%         set(gca,'YTick',[]); % Remove y axis ticks
%     end
title('Brain Graph vs. Creativity','FontSize',24);

F.fname=strcat(pre2, num2str(3));
F.wh=[3 2.5]*2;
print_fig(gcf,F)

% plot last figure
load(strcat(pre1,'CorrBrainNoiseSummary.mat'));
% cmap=zeros(3,3);
% ma = [1,0,1];
% cmap(1,:) = ma;
% cmap(2,:) = ma;
% map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

figure
%scatter(x,p(:,2), 500,'k.','jitter','on', 'jitterAmount', 0.3);
pv=p(:,2)
 [f,xi]=ksdensity(pv,'support',[0,1]); 
 hold on
 plot(xi,f,'.-','LineWidth',4);
 pv=sort(pv,'ascend');
 ord=0.01*ones(length(pv),1);
 for i=2:length(pv);
     if pv(i)-pv(i-1)<0.001
         ord(i)=ord(i-1)+0.4;
     end
 end
 plot(pv,ord,'.','MarkerSize',24);
xlim([0,0.15]);
ylim([-1 15]);
set(gca,'FontSize',24);
set(gca,'YTick',[]); % Remove y axis ticks
 hold off
xlabel('False Positive Rate','FontSize',24);
ylabel('Density Function','FontSize',24);
title('Brain Activity vs. Fake Movie','FontSize',24);

F.fname=strcat(pre2, 'CORR');
F.wh=[3 2.5]*2;
print_fig(gcf,F)
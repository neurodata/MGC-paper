function []=plot_realData(pre1,pre2)
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

%%
fpath = mfilename('fullpath');
findex=strfind(fpath,'\');
rootDir=fpath(1:findex(end-2));
p = genpath(rootDir);
gits=strfind(p,'.git');
colons=strfind(p,':');
for i=0:length(gits)-1
    endGit=find(colons>gits(end-i),1);
    p(colons(endGit-1):colons(endGit)-1)=[];
end
addpath(p);

if nargin<1
    pre1='..\..\Data\Results\'; % The folder to locate data
end
if nargin<2
    pre2='..\..\Figures\Fig'; % The folder to save figures
end
cmap=zeros(4,3);
gr = [0.5,0.5,0.5];
ma = [1,0,1];
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
filename=strcat(pre1,'CorrPermDistTestTypeBrainLMLxY.mat');
load(filename);
n=size(p1All,1);
% kmin=2;
% xa=kmin:n;
pp1=p1All(:,end);
filename=strcat(pre1,'CorrPermDistTestTypeBrainLMRxY.mat');
load(filename);
n=size(p1All,1);
pp2=p1All(:,end);

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

figure
hold on
ind1=1:n;
ind1=ind1(pp1>0.05);
ind1=[1 ind1 n+1];
[~,i]=max(ind1(2:end)-ind1(1:length(ind1)-1));
ind1=ind1(i)+1:ind1(i+1)-2;
ind1S=2:ind1(1);
ind1E=ind1(end):n;

ind2=1:n;
ind2=ind2(pp2>0.05);
ind2=[1 ind2 n+1];
[~,i]=max(ind2(2:end)-ind2(1:length(ind2)-1));
ind2=ind2(i)+1:ind2(i+1)-2;
ind2S=2:ind2(1);
ind2E=ind2(end):n;

plot(ind1S,pp1(ind1S),'-',ind2S,pp2(ind2S),'--',ind1E,pp1(ind1E),'-',ind2E,pp2(ind2E),'--','LineWidth',2);
h=legend('Left Brain','Right Brain','Location','North');
set(h,'FontSize',20);
legend boxoff
plot(ind1,pp1(ind1),'k-',ind2,pp2(ind2),'k--','LineWidth',4)
set(gca,'FontSize',16);
xlabel('Number of Neighbors for X','FontSize',20);
xlim([2 n]);
ylim([0 0.2]);
ylabel('P-Value','FontSize',20);

title('Brain Shape vs. Disorder','FontSize',24);

hold off
F.fname=strcat(pre2, num2str(1));
F.wh=[3 2.5]*2;
print_fig(gcf,F)

% Plot second figure
filename=strcat(pre1,'CorrPermDistTestTypeMigrainxCCI.mat');
load(filename);
figure
kmin=2;
imagesc(p1All(kmin:end,kmin:end)');
set(gca,'YDir','normal')
colormap(flipud(map2))
caxis([0.01 0.1])
set(gca,'FontSize',16);
%     if i==3
xlabel('Number of Neighbors for X','FontSize',20);
ylabel('Number of Neighbors for Y','FontSize',20);
colorbar
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
 [f,xi]=ksdensity(p(:,2),'support',[0,1]); 
 hold on
 plot(xi,f,'.-','LineWidth',2);
 plot(p(:,2),0.01*ones(size(p,1),1),'k.','MarkerSize',20);
xlim([0,0.15]);
ylim([-1 15]);
set(gca,'FontSize',16);
set(gca,'YTick',[]); % Remove y axis ticks
 hold off
xlabel('False Positive Rate','FontSize',20);
ylabel('Density Function','FontSize',20);
title('Brain Activity vs. Fake Movie','FontSize',24);

F.fname=strcat(pre2, 'CORR');
F.wh=[3 2.5]*2;
print_fig(gcf,F)
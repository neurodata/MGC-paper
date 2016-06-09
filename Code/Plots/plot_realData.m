function []=plot_realData(pre1,pre2)
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

if nargin<1
    pre1='../../Data/Results/'; % The folder to locate data
end
if nargin<2
    pre2='../../Draft/Figures/FigReal'; % The folder to save figures
end
cmap=zeros(3,3);
gr = [0.5,0.5,0.5];
ma = [1,0,1];
% cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = gr;
% cmap(3,:) = cy;
map1=cmap;
map2 = brewermap(128,'GnBu'); % brewmap

% Plot heatmap
total=4;
set(groot,'defaultAxesColorOrder',map1);
[filename, titleStr]=CorrRealDataName(1);
filename=strcat(pre1,filename);
load(filename);
n=size(p1All,1);
kmin=2;
xa=kmin:n;
pp1=p1All(xa,end);
[filename, titleStr]=CorrRealDataName(2);
filename=strcat(pre1,filename);
load(filename);
n=size(p1All,1);
pp2=p1All(xa,end);
figure

y=0.05;
H=area(xa,y*ones(n-kmin+1,1),'LineStyle',':','FaceColor', [.9 .9 .9]);
hold on
plot(xa,pp1,'k.-',xa,pp2,'.--','LineWidth',2)
set(gca,'FontSize',14);
xlabel('Number of Neighbors for X','FontSize',16);
xlim([2 n]);
ylim([0 0.2]);
ylabel('P-Value','FontSize',16);

title('Local Tests P-value for Brain Shape vs Disorder','FontSize',17);
legend('Significant P-value Area','P-value for Testing Left Brain','P-value for Testing Right Brain','Location','North');
legend boxoff

hold off
F.fname=strcat(pre2, num2str(1));
F.wh=[3 2.5]*2;
print_fig(gcf,F)

for i=3:total
    [filename, titleStr]=CorrRealDataName(i);
    filename=strcat(pre1,filename);
    load(filename);
    figure
    kmin=2;
    imagesc(p1All(kmin:end,kmin:end)');
    set(gca,'YDir','normal')
    colormap(flipud(map2))
    caxis([0.01 0.1])
    set(gca,'FontSize',14);
%     if i==3
        xlabel('Number of Neighbors for X','FontSize',16);
        ylabel('Number of Neighbors for Y','FontSize',16);
        colorbar
%     else
%         set(gca,'XTick',[]); % Remove x axis ticks
%         set(gca,'YTick',[]); % Remove y axis ticks
%     end
    
    % Figure title/labels
    titleStr = strcat('Local Tests P-value for ', titleStr);
    title(titleStr,'FontSize',17);
    
    F.fname=strcat(pre2, num2str(i));
    F.wh=[3 2.5]*2;
    print_fig(gcf,F)
end

load(strcat(pre1,'CorrBrainNoiseSummary.mat'));
cmap=zeros(3,3);
ma = [1,0,1];
cmap(1,:) = ma;
cmap(2,:) = ma;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

x=ones(size(p,1),1);
scatter(x,p(:,2), 500,'k.','jitter','on', 'jitterAmount', 0.3);
ylim([0,0.15]);
set(gca,'FontSize',16);
ax=gca;

set(gca,'XTick',[]); % Remove x axis ticks
ylabel('False Positive Rate','FontSize',16);
title('False Positive Rates for Brain vs Noise','FontSize',17);

F.fname=strcat(pre2, 'CORR');
F.wh=[3 2.5]*2;
print_fig(gcf,F)


function [str, title]=CorrRealDataName(i)
str='CorrPermDistTestType';
switch i
    %     case 1
    %         str=strcat(str,'BrainCxP.mat');
    %         title=' Connectome vs Personality';
    case 1
        str=strcat(str,'BrainLMLxY.mat');
        title=' Left Brain Shape vs Disorder';
    case 2
        str=strcat(str,'BrainLMRxY.mat');
        title=' Right Brain Shape vs Disorder';
    case 3
        str=strcat(str,'MigrainxCCI.mat');
        title=' Migrain vs CCI';
    case 4
        str=strcat(str,'M2gxCCI.mat');
        title=' M2g vs CCI';
end
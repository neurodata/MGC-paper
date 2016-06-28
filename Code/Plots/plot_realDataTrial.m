function []=plot_realDataTrial

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
map3 = brewermap(128,'PRGn'); % brewmap
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

filename=strcat(pre1,'CorrPermDistTestTypeBrainCxP.mat');
load(filename)
figNumber='CxP';
figure
kmin=2;
imagesc(p2All(kmin:end,kmin:end)');
set(gca,'YDir','normal')
colormap(flipud(map2))
caxis([0 0.1])
set(gca,'FontSize',24);
%     if i==3
[n1,n2]=size(p2All);
xlabel('# of Neighbors for X','FontSize',24);
ylabel('# of Neighbors for Y','FontSize',24);
set(gca,'XTick',[1,round(n1/2)-1,n1-1],'XTickLabel',[2,round(n1/2),n1]); % Remove x axis ticks
set(gca,'YTick',[1,round(n2/2)-1,n2-1],'YTickLabel',[2,round(n2/2),n2]); % Remove x axis ticks
h=colorbar('Ticks',[0,0.05,0.1]);%,'location','westoutside');
%     else
%         set(gca,'XTick',[]); % Remove x axis ticks
%         set(gca,'YTick',[]); % Remove y axis ticks
%     end
title('Brain vs. Personality','FontSize',24);

F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)
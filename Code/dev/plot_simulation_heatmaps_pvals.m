function []=plot_simulation_heatmaps_pvals(repPerm)
% Used to plot figure 1-8 used in tex. Run like

% total is usually 20.
% pre1 specifies the location to load data.
% pre2 specifies the location to save pictures.

%%% File path searching
if nargin<1
    repPerm=100;
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
pre2=strcat(rootDir,'Figures/Fig');% The folder to save figures

total=20;
map2 = brewermap(128,'GnBu'); % brewmap
map2=flipud(map2);
cticks=[0.01, 0.1, 0.5];
% set(groot,'defaultAxesColorOrder',map1);

figNumber='1DHeatP';
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
nn=60;
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=strcat(num2str(j),'.',{' '},CorrSimuTitle(j));
    
    [x, y]=CorrSampleGenerator(j,nn,1,1, 1);
    A=squareform(pdist(x));
    B=squareform(pdist(y));
    [~,~,ph]=MGCPermutationTest(A,B,repPerm,'mcor');
    ph=ph(2:end,2:end)';
    ph(ph<cticks(1))=cticks(1);
    ph(ph>cticks(end))=cticks(end);
%     % ind=[find(max(power2,[],1)>=thres,1) lim];
%     % lim=min(ind);
%     ind=find(numRange==nn);
%     if isempty(ind)
%         ind=1;
%     end
%     ph=powerMLocal(kmin:numRange(ind),kmin:numRange(ind),ind)';
    
    
    tt=find(sum(ph,2)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(tt:end,:)=repmat(ph(tt-1,:),numRange(ind)-tt,1);
    end
    tt=find(sum(ph,1)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(:,tt:end)=repmat(ph(:,tt-1),1,numRange(ind)-tt);
    end
    imagesc(log(ph));
    set(gca,'YDir','normal')
    colormap(map2)
%     caxis([cticks(1),cticks(3)])
    set(gca,'FontSize',14);
    set(gca,'XTick',[1,round(nn/2)-1,nn-1],'XTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
    set(gca,'YTick',[1,round(nn/2)-1,nn-1],'YTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
    title(titlechar,'FontSize',14, ...
        'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
    if j~=1
        set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
    end
    axis('square')
end
xlabel('# of X Neighbors','position',[-172 -12],'FontSize',24);
ylabel('# of Y Neighbors','position',[-430 156],'FontSize',24);
h=colorbar('Ticks',log(cticks),'TickLabels',cticks);
% h=colorbar('Ticks',[0,thres/2,thres]);
tstring=' of mcorr ';
h=suptitle(strcat('One-Dimensional Multiscale P-Value Maps'));
set(h,'FontSize',24,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)

figNumber='HDHeatP';
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=strcat(num2str(j),'.',{' '},CorrSimuTitle(j));
    kmin=2;thres=0.5;
    ind=find(max(powerMGCM,[],1)>=thres,1,'last');
    if isempty(ind)
        ind=1;
    end
    dd=dimRange(ind);
    
    [x, y]=CorrSampleGenerator(j,100,dd,1, 0);
    A=squareform(pdist(x));
    B=squareform(pdist(y));
    [~,~,ph]=MGCPermutationTest(A,B,repPerm,'mcor');
    ph=ph(2:end,2:end)';
    ph(ph<cticks(1))=cticks(1);
    ph(ph>cticks(end))=cticks(end);
%     ph=powerMLocal(kmin:n,kmin:n,ind)';
    tt=find(sum(ph,2)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(tt:end,:)=repmat(ph(tt-1,:),n-tt,1);
    end
    tt=find(sum(ph,1)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(:,tt:end)=repmat(ph(:,tt-1),1,n-tt);
    end
    imagesc(log(ph));
    set(gca,'YDir','normal')
    colormap(map2)
%     caxis([cticks(1),cticks(3)])
    set(gca,'XTick',[1,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n]); % Remove x axis ticks
    set(gca,'YTick',[1,round(n/2)-1,n-1],'YTickLabel',[2,round(n/2),n]); % Remove x axis ticks
    if j~=1
        set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
    end
    set(gca,'FontSize',14);
    title(titlechar,'FontSize',14);
    axis('square');
end
xlabel('# of X Neighbors','position',[-290 -20],'FontSize',24);
ylabel('# of Y Neighbors','position',[-720 260],'FontSize',24);
h=colorbar('Ticks',log(cticks),'TickLabels',cticks);
% h=colorbar('Ticks',[0,thres/2,thres]);
set(h,'FontSize',14);
h=suptitle(strcat('High-Dimensional Multiscale P-Value Maps'));
set(h,'FontSize',24,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)
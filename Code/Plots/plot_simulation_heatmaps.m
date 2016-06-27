function []=plot_simulation_heatmaps
% Used to plot figure 1-8 used in tex. Run like

% total is usually 20.
% pre1 specifies the location to load data.
% pre2 specifies the location to save pictures.

%%% File path searching
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
% set(groot,'defaultAxesColorOrder',map1);

figNumber='1DHeat';
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
        K=n;kmin=2;thres=0.8;
    %     ind=[find(max(power2,[],1)>=thres,1) lim];
    %     lim=min(ind);
    ind=find(numRange==50);
    if isempty(ind)
        ind=1;
    end
    if lim==21 && ind>1
        ind=ind-1;
    end
    ph=power2All(kmin:numRange(ind),kmin:numRange(ind),ind)';
    tt=find(sum(ph,2)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(tt:end,:)=repmat(ph(tt-1,:),numRange(ind)-tt,1);
    end
    tt=find(sum(ph,1)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(:,tt:end)=repmat(ph(:,tt-1),1,numRange(ind)-tt);
    end
    imagesc(ph);
    set(gca,'YDir','normal')
    colormap(map2)
    caxis([0 thres])
    set(gca,'FontSize',14);
    set(gca,'XTick',[1,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n]); % Remove x axis ticks
    set(gca,'YTick',[1,round(n/2)-1,n-1],'YTickLabel',[2,round(n/2),n]); % Remove x axis ticks
%     set(gca,'XTick',[]); % Remove x axis ticks
%     set(gca,'YTick',[]); % Remove y axis ticks
    title(titlechar,'FontSize',14);
    if j~=1
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    end
axis('square')
end
xlabel('# of Neighbors for X','position',[-290 -20],'FontSize',24);
ylabel('# of Neighbors for Y','position',[-720 300],'FontSize',24);
colorbar
h=colorbar('Ticks',[thres/2,thres]);%,'location','westoutside');
tstring=' of mcorr ';
h=suptitle(strcat('Multiscale Power Maps'));% for 1-Dimensional Simulations'));
set(h,'FontSize',24,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)

figNumber='HDHeat';
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    kmin=2;thres=0.5;
    ind=find(max(power2,[],1)>=thres,1,'last');
    if isempty(ind)
        ind=1;
    end
    if lim==21 && ind>1
        ind=ind-1;
    end
    ph=power2All(kmin:n,kmin:n,ind)';
 tt=find(sum(ph,2)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(tt:end,:)=repmat(ph(tt-1,:),n-tt,1);
    end
    tt=find(sum(ph,1)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(:,tt:end)=repmat(ph(:,tt-1),1,n-tt);
    end
%     if max(max(ph))>thres
%         ph=ph/max(max(ph))*thres; % in cas
%     end
    imagesc(ph);
    set(gca,'YDir','normal')
    colormap(map2)
    caxis([0 thres])
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
xlabel('# of Neighbors for X','position',[-290 -20],'FontSize',24);
ylabel('# of Neighbors for Y','position',[-720 300],'FontSize',24);
h=colorbar('Ticks',[thres/2,thres]);%,'location','westoutside');
set(h,'FontSize',14);
h=suptitle(strcat('Multiscale Power Maps'));% for High-Dimensional Simulations'));
set(h,'FontSize',24,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)
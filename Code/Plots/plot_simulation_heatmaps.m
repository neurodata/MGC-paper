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
    ind=[find(max(power2,[],1)>=thres,1) lim];
    lim=min(ind);
    ph=power2All(kmin:numRange(lim),kmin:numRange(lim),lim)';
 tt=find(sum(ph,2)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(tt:end,:)=repmat(ph(tt-1,:),numRange(lim)-tt,1);
    end
    tt=find(sum(ph,1)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(:,tt:end)=repmat(ph(:,tt-1),1,numRange(lim)-tt);
    end
    imagesc(ph);
    set(gca,'YDir','normal')
    colormap(map2)
    caxis([0 thres])
    set(gca,'FontSize',12);
%     set(gca,'XTick',[]); % Remove x axis ticks
%     set(gca,'YTick',[]); % Remove y axis ticks
    title(titlechar,'FontSize',14);
axis('square')
end
xlabel('# of Neighbors for X','position',[-290 -20],'FontSize',20);
ylabel('# of Neighbors for Y','position',[-720 300],'FontSize',20);
colorbar
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
    lim=ind;
    ph=power2All(kmin:n,kmin:n,lim)';
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
    if j~=16
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    end
    set(gca,'FontSize',12);
    title(titlechar,'FontSize',14);
    axis('square');
end
xlabel('# of Neighbors for X','position',[-290 -20],'FontSize',20);
ylabel('# of Neighbors for Y','position',[-720 300],'FontSize',20);
colorbar
h=suptitle(strcat('Multiscale Power Maps'));% for High-Dimensional Simulations'));
set(h,'FontSize',24,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)
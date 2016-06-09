function []=plot_outlierModel(pre1,pre2)
% Author: Cencheng Shen
% CorrVisualPlots()
% CorrVisualPlots(100,2)
% Used to plot figure 0 in the files
if nargin<1
    pre1='../../Data/Results/'; % The folder to locate data
end
if nargin<2
    pre2='../../Draft/Figures/Fig'; % The folder to save figures
end
cmap=zeros(3,3);
gr = [0.5,0.5,0.5];
ma = [1,0,1];
% cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = gr;
% cmap(3,:) = cy;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1)
map2 = brewermap(128,'GnBu'); % brewmap
%%% Outlier Model
seq=[0.3,0.5,0.7];
n=100;dim=1;

for i=1:length(seq)
    w=seq(i);
    figure;
    [x, y]=CorrSampleGenerator(0,n,dim,1,w);
    dep=1:ceil(n*w);
    ind=ceil(n*w)+1:n;
    plot(x(dep,1),y(dep,1),'k.',x(ind,1),y(ind,1),'.','MarkerSize',24);
    xlim([-5,10]);
    ylim([-10,5]);
    if i~=1
        set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
    else
        ylabel('Example Data','FontSize',24)%,'position',[-200 -0.2],'FontSize',20);
    end
    set(gca,'FontSize',16);
    title(strcat(num2str(w*100),'% outliers'),'FontSize',30);
    s1Pos=get(gca,'position');
    F.fname=strcat(pre2, 'OutlierVisual',num2str(i));
    F.wh=[3 2.5]*2;
    print_fig(gcf,F)
    
    
    figNumber=strcat('OutlierPower',num2str(i));
    filename=strcat(pre1,'CorrIndTestType0N100Dim1W',num2str(w),'.mat');
    load(filename)
    f2=figure;
    imagesc(power2All(2:end,2:end));
    set(gca,'YDir','normal')
    set(gca,'FontSize',16);
    colormap(map2)
    caxis([0 1])
    if i==3;
        h=colorbar('Ticks',[-1,0,1]);
        set(h,'FontSize',16,'location','eastoutside')       
        % Solution:
    end
    if i~=1
        set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
    else
        ylabel('Multiscale Power Map','FontSize',24)
    end
    s2Pos = get(gca,'position');
    s2Pos(3:4) = [s1Pos(3:4)];
    set(gca,'position',s2Pos);
    %title(strcat('Local Testing Powers at w=',num2str(w)),'FontSize',15);
    F.fname=strcat(pre2, figNumber);
    F.wh=[3 2.5]*2;
    print_fig(gcf,F)
end

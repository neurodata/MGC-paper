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

%% figure stuff

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
set(groot,'defaultAxesColorOrder',map1);


fnames={'CorrPermDistTestTypeBrainCxP.mat'; ...
    'CorrPermDistTestTypeBrainLMRxY.mat'; ...
    'CorrPermDistTestTypeMigrainxCCI.mat'};
xlabs={ 'Activity Neighbors'; ...
        'Shape Neighbors';...
        'Graph Neighbors'};
ylabs={ 'Personality Neighbors'; ...
        'Disease Neighbors';...
        'Creativity Neighbors'};
% tits= {'A. Brain Activity vs. Personality'; ...
%     'B. Brain Shape vs. Disorder';...
%     'C. Brain Graph vs. Creativity'};
tits={'A';'B';'C'};



%% loop maps
figure(1), clf, hold all
filename=strcat(pre1,fnames{2});
load(filename);
cmap2=flipud(map2);
fs=9;
cticks=[0.001, 0.01, 0.1, 0.5];
lw=2;

for i=1:3
    filename=strcat(pre1,fnames{i});
    load(filename);
    n=size(p2All,1);
    
    subplot(1,4,i)
    hold on
    imagesc(log(p2All'));
    set(gca,'YDir','normal')
    colormap(cmap2)
    set(gca,'FontSize',fs-1);
    %     if i==3
    [n1,n2]=size(p2All);
    xlabel(xlabs{i},'FontSize',fs, ...
       'Units', 'normalized','Position', [-0.01, -0.14], 'HorizontalAlignment', 'left')
    ylabel(ylabs{i},'FontSize',fs, ...
        'Units', 'normalized','Position', [-0.16 0], 'HorizontalAlignment', 'left')
    title(tits{i},'FontSize',fs+1, ...
        'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')
    
    
    %[~,indP]=MGCScaleVerify(p2All',rep);
    indP=ind2;
    if indP(end)==size(p2All,1)*size(p2All,2)
        indP2=indP(1:end-1);
    else
        indP2=indP;
    end
    [J,I]=ind2sub(size(p2All'),indP2);
    Ymin=min(I);
    Ymax=max(I);
    Xmin=min(J);
    Xmax=max(J);
    
    plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
    plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
    plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
    plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
    if indP(end)==size(p2All,1)*size(p2All,2)
        plot(size(p2All,1),size(p2All,2),'g.','markerSize',16)
    end
    xticks=[5,round(n1/2)-1,n1-1];
    if i==1,  xticks(1)=3; end
    set(gca,'XTick',xticks,'XTickLabel',[2,round(n1/2),n1]); % Remove x axis ticks
    set(gca,'YTick',[3,round(n2/2)-1,n2-1],'YTickLabel',[2,round(n2/2),n2]); % Remove x axis ticks
    xlim([2,n1]);
    ylim([2,n2]);
    axis('square');
    hold off
    if i==1
        hc=colorbar('Ticks',log(cticks),'TickLabels',cticks,'location','westoutside','FontSize',fs+1);
        title(hc,'p-value')
    end
end

% plot last figure
load(strcat(pre1,'CorrBrainNoiseSummary.mat'));
% cmap=zeros(3,3);
% ma = [1,0,1];
% cmap(1,:) = ma;
% cmap(2,:) = ma;
% map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

subplot(1,4,4)
%scatter(x,p(:,2), 500,'k.','jitter','on', 'jitterAmount', 0.3);
pv=p(:,2);
[f,xi]=ksdensity(pv);
hold on
plot(xi,f,'.-','LineWidth',lw);
pv=sort(pv,'ascend');
ord=0.01*ones(length(pv),1);
for i=2:length(pv);
    if pv(i)-pv(i-1)<0.001
        ord(i)=ord(i-1)+0.4;
    end
end
plot(pv,ord,'.','MarkerSize',8);
% <<<<<<< HEAD
xlim([0,0.15]);
ylim([-1 max(f)+1]);
% =======
% xlim([-0.05,0.15]);
% ylim([-1 15]);
% >>>>>>> 7132b2753b089edc2ca608140c119b901e31b17a
set(gca,'FontSize',fs-2);
set(gca,'YTick',[]); % Remove y axis ticks
axis('square');
hold off
xlabel('False Positive Rate','FontSize',fs,...
    'Units', 'normalized','Position', [-0.03, -0.12], 'HorizontalAlignment', 'left')
ylabel('Density','FontSize',fs, ...
    'Units', 'normalized','Position', [-0.05 0], 'HorizontalAlignment', 'left')
% title('D. Brain Activity vs. Fake Movie','FontSize',fs+2, ...
%     'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')
title('D','FontSize',fs+2, ...
    'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')

% F.fname=strcat(pre2, 'CORR');
% F.wh=[3 2.5]*2;
% print_fig(gcf,F)
%colorbar()

h=suptitle(strcat('Brain vs Mental Properties'));% for 1-Dimensional Simulations'));
set(h,'FontSize',15,'FontWeight','normal');
        
F.fname=pre2; %strcat(pre2, num2str(i));
F.wh=[8.2 2.2];
print_fig(gcf,F)

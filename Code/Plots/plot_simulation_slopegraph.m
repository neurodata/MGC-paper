function []=plot_simulation_slopegraph
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
fontSize=20;
%% Set colors
map1=zeros(7,3);
gr = [0.5,0.5,0.5];
cmap(1,:) = gr;
cmap(2,:) = gr;
cmap(3,:) = gr;
cmap(7,:) = [0 0 0];
dcorr = cmap(1,:);
mcorr = cmap(2,:);
mante = cmap(3,:);
HHG   = [0.5,0.5,0.5];
map1(1,:)=mante; map1(5,:)=mante; % The color for MGC{dcorr} and global dcorr.
map1(2,:)=dcorr; map1(6,:)=dcorr; % The color for MGC{mcorr} and global mcorr.
map1(3,:)=mcorr; map1(7,:)=mcorr; % The color for MGC{Mantel} and global Mantel.
map1(4,:)=mcorr; map1(8,:)=mcorr; % The color for MGC{Mantel} and global Mantel.
set(groot,'defaultAxesColorOrder',map1);

figure('units','normalized','position',[0 0 1 1])
s=2;t=4;
ls{3}='-';
ls{2}='--';
ls{1}=':';
ls{4}='-';

AUC=zeros(8,20);
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    AUC(1,j)=mean(powerP);
    AUC(2,j)=mean(powerMGCP);
    AUC(3,j)=mean(powerD);
    AUC(4,j)=mean(powerMGCD);
    AUC(5,j)=mean(powerM);
    AUC(6,j)=mean(powerMGCM);
    AUC(7,j)=mean(powerM);
    AUC(8,j)=mean(powerMGC);
end

x=1:2;
str={'A. Mantel';'B. Dcorr';'C. Mcorr';'D. Mcorr'};
for j=1:4
    subplot(s,t,j);
    hold on
    for i=6:19
        plot(x,AUC(2*j-1:2*j,i),ls{j},'LineWidth',1,'Color',map1(j,:));
    end
    plot(x,[mean(AUC(2*j-1,6:19)) mean(AUC(2*j,6:19))],ls{j},'LineWidth',10,'Color',map1(j,:));
    ylim([0,1]);
    set(gca,'XTick',[]); % Remove y axis ticks
    set(gca,'YTickLabel',[0,1],'YTick',0:1,'FontSize',fontSize);
    if j==1
        ylabel('1-Dimensional','FontSize',fontSize);
    else
        set(gca,'YTick',[]); % Remove y axis ticks
    end
    title(str(j,:),'FontSize',fontSize, ...
        'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
    hold off
end


AUC=zeros(8,20);
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    AUC(1,j)=mean(powerP);
    AUC(2,j)=mean(powerMGCP);
    AUC(3,j)=mean(powerD);
    AUC(4,j)=mean(powerMGCD);
    AUC(5,j)=mean(powerM);
    AUC(6,j)=mean(powerMGCM);
    AUC(7,j)=mean(powerM);
    AUC(8,j)=mean(powerMGC);
end

strL={'Global';'OMGC'};
for j=1:4
    subplot(s,t,4+j);
    hold on
    for i=6:19
        plot(x,AUC(2*j-1:2*j,i),ls{j},'LineWidth',1,'Color',map1(j,:));
    end
    plot(x,[mean(AUC(2*j-1,6:19)) mean(AUC(2*j,6:19))],ls{j},'LineWidth',10,'Color',map1(j,:));
    ylim([0,1]);
    if j==4
        strL={'Global';'SMGC'};
    end
    set(gca,'XTickLabel',strL,'XTick',[1.1,1.9],'YTickLabel',[0,1],'YTick',0:1,'FontSize',fontSize,'TickLength',[0 1]);
    if j==1
        ylabel('High-Dimensional','FontSize',fontSize);
    else
        set(gca,'YTick',[]); % Remove y axis ticks
    end
    hold off
end
F.fname=strcat(pre2,'Slope');
F.wh=[6 2.4]*2;
print_fig(gcf,F)

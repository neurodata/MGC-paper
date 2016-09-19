function []=plot_simulation_powerCompare
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
fontSize=16;
%% Set colors
map1=zeros(7,3);
gr = [0.5,0.5,0.5];
bl=[0,0,0];
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
s=1;t=4;
ls{3}='-';
ls{2}='--';
ls{1}=':';
ls{4}='-';

AUC=zeros(8,20);
% %load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    AUC(1,j)=mean(powerMGCP)-mean(powerP);
    AUC(2,j)=mean(powerMGCD)-mean(powerD);
    AUC(3,j)=mean(powerMGCM)-mean(powerM);
    AUC(4,j)=mean(powerMGC)-mean(powerM);
end
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    AUC(5,j)=mean(powerMGCP)-mean(powerP);
    AUC(6,j)=mean(powerMGCD)-mean(powerD);
    AUC(7,j)=mean(powerMGCM)-mean(powerM);
    AUC(8,j)=mean(powerMGC)-mean(powerM);
end
% 
str={'A. Mantel';'B. Dcorr';'C. Mcorr';'D. Mcorr'};
lw=3;
for j=1:4
    subplot(s,t,j);
    [f,xi]=ksdensity(AUC(j,6:19));
    hold on
    plot(xi,f,'.-','LineWidth',lw,'Color',gr);
    [f,xi]=ksdensity(AUC(4+j,6:19));
    plot(xi,f,'.-','LineWidth',lw,'Color',bl);
    if j==1
    legend('1-D','H-D','Position','NorthEast');
    end
    
    pv=sort(AUC(j,6:19),'ascend');
    ord=0.01*ones(length(pv),1);
    for i=2:length(pv);
        if pv(i)-pv(i-1)<0.005
            ord(i)=ord(i-1)+0.1;
        end
    end
    plot(pv-0.01,ord,'.','MarkerSize',15,'Color',gr);
    
    pv=sort(AUC(4+j,6:19),'ascend');
    ord=0.01*ones(length(pv),1);
    for i=2:length(pv);
        if pv(i)-pv(i-1)<0.005
            ord(i)=ord(i-1)+0.1;
        end
    end
    plot(pv+0.01,ord,'.','MarkerSize',20,'Color',bl);
    
    axis('square');
    hold off
    xlim([-0.1,1]);
    ylim([0,3]); 
    set(gca,'FontSize',fontSize);
    if j==1
        xlabel('Power Difference','FontSize',fontSize,...
            'Units', 'normalized','Position', [0.55, -0.2]);%, 'HorizontalAlignment', 'left')
        ylabel('Density','FontSize',fontSize, ...
           'Units', 'normalized','Position', [-0.2 0.23])%, 'HorizontalAlignment', 'left')
        set(gca,'XTickLabel',[0,1],'XTick',[0,1],'FontSize',fontSize);
        set(gca,'YTickLabel',[0,3],'YTick',[0,3],'FontSize',fontSize);
    else
        set(gca,'YTick',[]); % Remove y axis ticks
        set(gca,'XTick',[]); % Remove y axis ticks
    end
    title(str(j,:),'FontSize',fontSize, ...
       'Units', 'normalized','Position', [0.3 1.05])%, 'HorizontalAlignment', 'left')
    hold off
end

F.png=1;
F.fname=strcat(pre2,'Slope');
F.wh=[5 1.5]*2;
print_fig(gcf,F)

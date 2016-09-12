function []=plot_simulation_powers(select)
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

%%
if nargin<1
    select=1;
end
total=20;

lowd=0;
highd=1;

%% Set colors
map2 = brewermap(128,'PiYG'); % brewmap
%loca=map2(110,:);
%glob=map2(18,:);
loca=[0,1,0];
glob= [1,0,1];
HHG   = [0.5,0.5,0.5];

%set(groot,'defaultAxesColorOrder',map1);

ls{1}='-';
ls{2}='--';
ls{3}=':';
%ls{4}='.:';
ls{4}='--';

close all

%% figure1-4
figNumber='1DPower';
if select~=1
    figNumber='1DPowerAll';
end
figure('units','normalized','position',[0 0 1 1])
%set(groot,'defaultAxesColorOrder',map1);
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    %
    hold on
    if select==1
        h2=plot(numRange,powerM,ls{1},'LineWidth',3,'Color',glob);
        h1=plot(numRange,powerMGCM,ls{1},'LineWidth',3,'Color',loca);
        h4=plot(numRange,powerMGC,ls{1},'LineWidth',3,'Color','cyan');
        h3=plot(numRange,powerHHG,ls{4},'LineWidth',3,'Color',HHG);
    else
        h6=plot(numRange,powerP,ls{3},'LineWidth',3,'Color',glob);
        h5=plot(numRange,powerD,ls{2},'LineWidth',3,'Color',glob);
        h4=plot(numRange,powerM,ls{1},'LineWidth',3,'Color',glob);
        h3=plot(numRange,powerMGCP,ls{3},'LineWidth',3,'Color',loca);
        h2=plot(numRange,powerMGCD,ls{2},'LineWidth',3,'Color',loca);
        h1=plot(numRange,powerMGCM,ls{1},'LineWidth',3,'Color',loca);
        h8=plot(numRange,powerMGC,ls{1},'LineWidth',3,'Color','cyan');
        h7=plot(numRange,powerHHG,ls{4},'LineWidth',3,'Color',HHG);
    end
    hold off
    xlim([numRange(1) numRange(end)]);
    ylim([0 1]);
    if j~=1 % Remove x&y axis ticks except type 16, which is at the left bottom
        set(gca,'YTick',[]); % Remove y axis ticks
        set(gca,'XTick',[]); % Remove x axis ticks
    else
        set(gca,'XTick',[numRange(1),numRange(end)]); % Remove x axis ticks
        set(gca,'YTick',[0,0.5,1]); % Remove x axis ticks
    end
    set(gca,'FontSize',14);
    %    set(gca,'XTickLabel','FontSize',16);
    title(titlechar,'FontSize',14, ...
        'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
    axis('square');
end
xlabel('Sample Size','position',[-270 -0.2],'FontSize',24);
ylabel('Power','position',[-687 2.7],'FontSize',24);
h=suptitle('Testing Power for 20 Simulated 1-Dimensional Settings');
set(h,'FontSize',24,'FontWeight','normal');
lgdPosition = [0.03, 0.85, .05, .05]; %Legend Position
if select==1;
    % h=legend([h1 h2 h3],'MGC','Mcorr','HHG','Location',lgdPosition);
    h=legend([h1 h4 h2 h3],'Oracle MGC','Sample MGC','Mcorr', 'HHG','Location',lgdPosition);
else
    h=legend([h1 h2 h3 h8 h4 h5 h6 h7],'MGC_{M}','MGC_{D}','MGC_{P}','Sample MGC','Mcorr','Dcorr','Mantel','HHG','Location',lgdPosition);
end
legend boxoff
set(h,'FontSize',14);
%
F.fname=[strcat(pre2, figNumber)]; %, '_', num2str(cm)];
F.wh=[8 5]*2;
print_fig(gcf,F)

%% Plot 5-8
% if select==1;
% map1(1,:)=mcorr; map1(3,:)=mcorr; % The color for MGC{dcorr} and global dcorr.
% map1(2,:)=HHG; % The color for HHG
% end

figNumber='HDPower';
if select~=1
    figNumber='HDPowerAll';
end
figure('units','normalized','position',[0 0 1 1])
%set(groot,'defaultAxesColorOrder',map1)
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    numRange=dimRange;
    subplot(s,t,j)
    titlechar=[CorrSimuTitle(j)]; %, ' d=', num2str(max(dimRange))];
    hold on
    if select==1
        h2=plot(numRange,powerM,ls{1},'LineWidth',3,'Color',glob);
        h1=plot(numRange,powerMGCM,ls{1},'LineWidth',3,'Color',loca);
        h4=plot(numRange,powerMGC,ls{1},'LineWidth',3,'Color','cyan');
        h3=plot(numRange,powerHHG,ls{4},'LineWidth',3,'Color',HHG);
    else
        h6=plot(numRange,powerP,ls{3},'LineWidth',3,'Color',glob);
        h5=plot(numRange,powerD,ls{2},'LineWidth',3,'Color',glob);
        h4=plot(numRange,powerM,ls{1},'LineWidth',3,'Color',glob);
        h3=plot(numRange,powerMGCP,ls{3},'LineWidth',3,'Color',loca);
        h2=plot(numRange,powerMGCD,ls{2},'LineWidth',3,'Color',loca);
        h1=plot(numRange,powerMGCM,ls{1},'LineWidth',3,'Color',loca);
        h8=plot(numRange,powerMGC,ls{1},'LineWidth',3,'Color','cyan');
        h7=plot(numRange,powerHHG,ls{4},'LineWidth',3,'Color',HHG);
    end
    hold off
    xlim([numRange(1) numRange(end)]);
    ylim([0 1]);
    if j~=1 % Remove x&y axis ticks except type 16, which is at the left bottom
        set(gca,'YTick',[]); % Remove y axis ticks
    else
        set(gca,'YTick',[0,0.5,1]); % Remove x axis ticks
    end
    set(gca,'XTick',[numRange(1),numRange(end)]); % Remove x axis ticks
    set(gca,'FontSize',14);
    title(titlechar,'FontSize',14, ...
        'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
    %set(gca,'XTickLabel','FontSize',16);
    axis('square');
end
xlabel('Dimension','position',[-290 -0.2],'FontSize',24);
ylabel('Power','position',[-720 2.7],'FontSize',24);
h=suptitle('Testing Power for 20 Simulated High-Dimensional Settings');
set(h,'FontSize',24,'FontWeight','normal');
lgdPosition = [0.03, 0.85, .05, .05]; %Legend Position
if select==1;
    %h=legend([h1 h2 h3],'MGC','Mcorr','HHG','Location',lgdPosition);
    h=legend([h1 h4 h2 h3],'Oracle MGC','Sample MGC','Mcorr', 'HHG','Location',lgdPosition);
else
    h=legend([h1 h2 h3 h8 h4 h5 h6 h7],'MGC_{M}','MGC_{D}','MGC_{P}','Sample MGC','Mcorr','Dcorr','Mantel','HHG','Location',lgdPosition);
end
set(h,'FontSize',14);
legend boxoff
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)
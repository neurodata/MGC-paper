function []=plot_simulation_powers2(select)
% Used to plot figure 1-8 used in tex. Run like

% total is usually 20.
% pre1 specifies the location to load data.
% pre2 specifies the location to save pictures.

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
pre2=strcat(rootDir,'Figures/Fig');% The folder to save figures

%%
if nargin<1
    select=1;
end
total=20;

%% Set colors
map3 = brewermap(128,'PiYG'); % brewmap
mgc='green';
loca='cyan';
glob= [0.5,0.5,0.5];
HHG   = 'magenta';
mcorr='cyan';
mantel='cyan';
dcorr='cyan';
pcorr=[0.8,0.8,0.8];
mic   = [0.3,0.3,0.3];
kendall=[0.8,0.8,0.8];
spearman=[0.8,0.8,0.8];
hsic='blue';

lw=3;
ls{1}='-';
ls{2}='--';
ls{3}='-.';
ls{4}=':';

%% figure1-4
figNumber='1DPowerMantel';
if select~=1
    figNumber='1DPowerAll';
end
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=strcat(num2str(j),'.',{' '},CorrSimuTitle(j));
    %
    hold on
%     if select==1
        h7=plot(numRange,powerP-powerMGC3,ls{3},'LineWidth',3,'Color',mantel);
%         h6=plot(numRange,powerD-powerMGC,ls{2},'LineWidth',3,'Color',dcorr);
        h5=plot(numRange,powerM-powerMGC,ls{1},'LineWidth',3,'Color',mgc);
        h0=plot(numRange,zeros(length(numRange)),ls{4},'LineWidth',1,'Color',glob);
        %h4=plot(numRange,powerMGCP,ls{3},'LineWidth',3,'Color',loca);
        %h3=plot(numRange,powerMGCD,ls{2},'LineWidth',3,'Color',loca);
        %h2=plot(numRange,powerMGCM-powerMGC,ls{1},'LineWidth',3,'Color',loca);
%         h1=plot(numRange,powerMGC-powerMGC,ls{1},'LineWidth',3,'Color',MGC);
%         h8=plot(numRange,powerHHG-powerMGC,ls{4},'LineWidth',3,'Color',HHG);
%         h9=plot(numRange,powerHSIC-powerMGC,ls{3},'LineWidth',3,'Color',hsic);
%         h10=plot(numRange,powerCorr-powerMGC,ls{2},'LineWidth',3,'Color',pcorr);
%         h11=plot(numRange,powerMIC-powerMGC,ls{1},'LineWidth',3,'Color',mic);
%     else
%         h7=plot(numRange,powerP,ls{3},'LineWidth',3,'Color',glob);
%         h6=plot(numRange,powerD,ls{2},'LineWidth',3,'Color',glob);
%         h5=plot(numRange,powerM,ls{1},'LineWidth',3,'Color',glob);
%         h4=plot(numRange,powerMGCP,ls{3},'LineWidth',3,'Color',loca);
%         h3=plot(numRange,powerMGCD,ls{2},'LineWidth',3,'Color',loca);
%         h2=plot(numRange,powerMGCM,ls{1},'LineWidth',3,'Color',loca);
%         h1=plot(numRange,powerMGC,ls{1},'LineWidth',3,'Color',MGC);
%         h8=plot(numRange,powerHHG,ls{4},'LineWidth',3,'Color',HHG);
%         h9=plot(numRange,powerHSIC,ls{3},'LineWidth',3,'Color',hsic);
%         h10=plot(numRange,powerCorr,ls{2},'LineWidth',3,'Color',pcorr);
%         h11=plot(numRange,powerMIC,ls{1},'LineWidth',3,'Color',mic);
%     end
    hold off
    xlim([numRange(1) numRange(end)]);
    ylim([-1 1]);
    if j~=1 % Remove x&y axis ticks except type 16, which is at the left bottom
        set(gca,'YTick',[]); % Remove y axis ticks
        set(gca,'XTick',[]); % Remove x axis ticks
    else
        set(gca,'XTick',[numRange(1),numRange(end)]); % Remove x axis ticks
        set(gca,'YTick',[-1,-0.5,0,0.5,1]); % Remove x axis ticks
    end
    set(gca,'FontSize',14);
    title(titlechar,'FontSize',14, ...
        'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
    axis('square');
end
% xlabel('Sample Size','position',[-270 -0.2],'FontSize',24);
% ylabel('Power','position',[-687 2.7],'FontSize',24);
xlabel('Sample Size','position',[-270 -1.6],'FontSize',24);
ylabel('Power Difference','position',[-687 4.5],'FontSize',24);
h=suptitle('Testing Power Difference for 1D Settings');
set(h,'FontSize',24,'FontWeight','normal');
lgdPosition = [0.05, 0.92, .05, .05]; %Legend Position
% if select==1;
    h=legend([h5 h7],'Mcorr - MGC (Mcorr)','Mantel - MGC(Mantel)','Location',lgdPosition);
% else
%     h=legend([h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11],'Sample MGC','MGC_{M}','MGC_{D}','MGC_{P}','Mcorr','Dcorr','Mantel','HHG','HSIC','Pearson','MIC','Location',lgdPosition);
% end
legend boxoff
set(h,'FontSize',14);
% %
F.fname=[strcat(pre2, figNumber)];
F.wh=[8 5]*2;
print_fig(gcf,F)

%% Plot 5-8
figNumber='HDPowerMantel';
if select~=1
    figNumber='HDPowerAll';
end
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    numRange=dimRange;
    subplot(s,t,j)
    titlechar=strcat(num2str(j),'.',{' '},CorrSimuTitle(j));
    hold on
%     if select==1
        h7=plot(numRange,powerP-powerMGC3,ls{3},'LineWidth',3,'Color',mantel);
%         h6=plot(numRange,powerD-powerMGC,ls{2},'LineWidth',3,'Color',dcorr);
        h5=plot(numRange,powerM-powerMGC,ls{1},'LineWidth',3,'Color',mgc);
        h0=plot(numRange,zeros(length(numRange)),ls{4},'LineWidth',1,'Color',glob);
        %h4=plot(numRange,powerMGCP,ls{3},'LineWidth',3,'Color',loca);
        %h3=plot(numRange,powerMGCD,ls{2},'LineWidth',3,'Color',loca);
        %h2=plot(numRange,powerMGCM-powerMGC,ls{1},'LineWidth',3,'Color',loca);
%         h1=plot(numRange,powerMGC-powerMGC,ls{1},'LineWidth',3,'Color',MGC);
%         h8=plot(numRange,powerHHG-powerMGC,ls{4},'LineWidth',3,'Color',HHG);
%         h9=plot(numRange,powerHSIC-powerMGC,ls{4},'LineWidth',3,'Color',hsic);
%         h10=plot(numRange,powerCorr-powerMGC,ls{4},'LineWidth',3,'Color',pcorr);
%         h11=plot(numRange,powerCCA-powerMGC,ls{3},'LineWidth',3,'Color',pcorr);
%     else
%         h7=plot(numRange,powerP,ls{3},'LineWidth',3,'Color',glob);
%         h6=plot(numRange,powerD,ls{2},'LineWidth',3,'Color',glob);
%         h5=plot(numRange,powerM,ls{1},'LineWidth',3,'Color',glob);
%         h4=plot(numRange,powerMGCP,ls{3},'LineWidth',3,'Color',loca);
%         h3=plot(numRange,powerMGCD,ls{2},'LineWidth',3,'Color',loca);
%         h2=plot(numRange,powerMGCM,ls{1},'LineWidth',3,'Color',loca);
%         h1=plot(numRange,powerMGC,ls{1},'LineWidth',3,'Color',MGC);
%         h8=plot(numRange,powerHHG,ls{4},'LineWidth',3,'Color',HHG);
%         h9=plot(numRange,powerHSIC,ls{4},'LineWidth',3,'Color',hsic);
%         h10=plot(numRange,powerCorr,ls{4},'LineWidth',3,'Color',pcorr);
%         h11=plot(numRange,powerCCA,ls{3},'LineWidth',3,'Color',pcorr);
%     end
    hold off
    xlim([numRange(1) numRange(end)]);
    ylim([-1 1]);
    if j~=1 % Remove x&y axis ticks except type 16, which is at the left bottom
        set(gca,'YTick',[]); % Remove y axis ticks
    else
        set(gca,'YTick',[-1,-0.5,0,0.5,1]); % Remove x axis ticks
    end
    set(gca,'XTick',[numRange(1),numRange(end)]); % Remove x axis ticks
    set(gca,'FontSize',14);
    title(titlechar,'FontSize',14, ...
        'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
    axis('square');
end
%xlabel('Dimension','position',[-290 -0.2],'FontSize',24);
%ylabel('Power','position',[-720 2.7],'FontSize',24);
xlabel('Dimension','position',[-290 -1.6],'FontSize',24);
ylabel('Power Difference','position',[-720 4.5],'FontSize',24);
h=suptitle('Testing Power Difference for HD Settings');
set(h,'FontSize',24,'FontWeight','normal');
lgdPosition = [0.05, 0.92, .05, .05]; %Legend Position
% if select==1;
h=legend([h5 h7],'Mcorr - MGC (Mcorr)','Mantel - MGC(Mantel)','Location',lgdPosition);
%     h=legend([h1 h6 h5 h7 h8 h9 h10 h11],'MGC','Dcorr','Mcorr','Mantel','HHG','HSIC','RV','CCA','Location',lgdPosition);
% else
%     h=legend([h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11],'Sample MGC','MGC_{M}','MGC_{D}','MGC_{P}','Mcorr','Dcorr','Mantel','HHG','HSIC','Pearson','CCA','Location',lgdPosition);
% end
set(h,'FontSize',14);
legend boxoff
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)
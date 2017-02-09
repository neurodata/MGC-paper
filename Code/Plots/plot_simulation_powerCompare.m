function []=plot_simulation_powerCompare(select)
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

total=20;
if nargin<1
    select=1;
end
fontSize=20;
%% Set colors
gr = [0.5,0.5,0.5];
bl=[0,0,0];

figure('units','normalized','position',[0 0 1 1])
s=1;t=4;
ls{3}='-';
ls{2}='--';
ls{1}=':';
%ls{4}='-';

AUC=zeros(8,total);
% %load data
cons=0;
if select==0
    for j=1:total
        filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
        load(filename)
        AUC(1,j+cons)=mean(powerMGCM)-mean(powerP);
        AUC(2,j+cons)=mean(powerMGCM)-mean(powerD);
        AUC(3,j+cons)=mean(powerMGCM)-mean(powerM);
        AUC(4,j+cons)=mean(powerMGCM)-mean(powerHHG);
        AUC(5,j+cons)=mean(powerMGC)-mean(powerP);
        AUC(6,j+cons)=mean(powerMGC)-mean(powerD);
        AUC(7,j+cons)=mean(powerMGC)-mean(powerM);
        AUC(8,j+cons)=mean(powerMGC)-mean(powerHHG);
    end
else
    for j=1:total
        filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
        load(filename)
        AUC(1,j+cons)=mean(powerMGCM)-mean(powerP);
        AUC(2,j+cons)=mean(powerMGCM)-mean(powerD);
        AUC(3,j+cons)=mean(powerMGCM)-mean(powerM);
        AUC(4,j+cons)=mean(powerMGCM)-mean(powerHHG);
        AUC(5,j+cons)=mean(powerMGC)-mean(powerP);
        AUC(6,j+cons)=mean(powerMGC)-mean(powerD);
        AUC(7,j+cons)=mean(powerMGC)-mean(powerM);
        AUC(8,j+cons)=mean(powerMGC)-mean(powerHHG);
    end
end
%
str={'Mantel';'Dcorr';'Mcorr';'HHG'};
%lw=3;
%ss=[1,3,2,5,4,6,7,13,19,16,17,10,9,14,8,11,12,15,18,20];
for j=1:length(str)
    subplot(s,t,j);
    hold on
    pv1=AUC(j,:);
    %pv1=pv1(ss);
    ord1=0.05*(1:total);
%     plot(pv2,ord1,'*','MarkerSize',8,'Color',gr);
    
    pv2=AUC(j+4,:);
    %pv2=pv2(ss);
    plot(pv2,ord1,'.','MarkerSize',15,'Color',bl);
        plot(pv2,ord1,'*','MarkerSize',8,'Color',bl);
%     if j==1
%         h=legend('Oracle','Sample');
%         set(h,'FontSize',fontSize, 'Units', 'normalized','Position',[0.06,0.75,0,0]);
%         legend boxoff
%     end
%     plot(pv1,ord1,'.','MarkerSize',15,'Color',gr);
    
    plot(zeros(21,1),0:0.05:1,'--','Color',gr,'LineWidth',2)
    mm1=mean(pv1);
    mm2=mean(pv2);
    pos=0;
%     plot(mm1,pos,'.','MarkerSize',30,'Color',gr);
    plot(mm2,pos,'*','MarkerSize',13,'Color',bl);
    plot(mm2,pos,'.','MarkerSize',24,'Color',bl);
    
%     a=text(mm1,pos+0.02,num2str(round(mm1*100)/100),'VerticalAlignment','bottom','HorizontalAlignment','left','Color',gr,'Interpreter','latex');
%     set(a,'FontSize',fontSize);
%     b=text(mm2,pos-0.03,num2str(round(mm2*100)/100),'VerticalAlignment','top','HorizontalAlignment','left','Color',bl,'Interpreter','latex');
%     set(b,'FontSize',fontSize);
    
    axis('square');
    hold off
    xlim([-1,1]);
    ylim([0,1]);
    set(gca,'FontSize',fontSize);
    set(gca,'XTick',[]); % Remove y axis ticks
    %set(gca,'XTickLabel',[0],'XTick',[0],'FontSize',fontSize);
    labels={-1,[],[],[],0,[],[],[],1};
    set(gca,'XTickLabel',labels,'XTick',-1:0.25:1,'FontSize',fontSize);
    yt=[1,5,10,15,20];
    set(gca,'YTickLabel',yt,'YTick',0.05*yt,'FontSize',fontSize);
        
    if j==1
        xlabel('Power Difference','FontSize',fontSize,...
            'Units', 'normalized','Position', [0.4, -0.18]);%, 'HorizontalAlignment', 'left')
        ylabel('Simulation Type','FontSize',fontSize, ...
            'Units', 'normalized','Position', [-0.18 0.38])%, 'HorizontalAlignment', 'left')
        %set(gca,'XTickLabel',[0,1],'XTick',[0,1],'FontSize',fontSize);
    else
        set(gca,'YTick',[]); % Remove y axis ticks
        %set(gca,'XTick',[]); % Remove y axis ticks
    end
    %set(gca,'YTick',[]); % Remove y axis ticks
    
%     switch j
%         case 1
%             tit=strcat('A. MGC - ',{' '}, str(j,:));
%         case 2
%             tit=strcat('B. MGC - ',{' '}, str(j,:));
%         case 3
%             tit=strcat('C. MGC - ',{' '}, str(j,:));
%         case 4
%             tit=strcat('D. MGC - ',{' '}, str(j,:));
%     end
    switch j
        case 1
            tit=strcat('A.');
        case 2
            tit=strcat('B.');
        case 3
            tit=strcat('C.');
        case 4
            tit=strcat('D.');
    end

    txt1=strcat(tit,' Power(MGC)');
    txt2=strcat({'     '},'- Power(',str(j,:),')');
    txt2=txt2{1};
    title({txt1,txt2},'FontSize',fontSize, ...
        'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
    hold off
end

if select~=1
    figNumber='1DPowerMGCM';
else
    figNumber='HDPowerMGCM';
end
F.fname=strcat(pre2,figNumber);
F.wh=[8 2.2]*2;
print_fig(gcf,F)

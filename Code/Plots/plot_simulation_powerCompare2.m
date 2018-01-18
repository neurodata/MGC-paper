function []=plot_simulation_powerCompare2
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
fontSize=18;
%% Set colors
% loca=[0,1,0];
% glob= [0.5,0.5,0.5];
HHG   = 'magenta';
mcorr='cyan';
mantel='cyan';
dcorr='cyan';
map3 = brewermap(128,'PiYG'); % brewmap
MGC=map3(100,:);
pcorr=[0.6,0.6,0.6];
mic   = [0.3,0.3,0.3];
kendall=[0.6,0.6,0.6];
spearman=[0.6,0.6,0.6];
hsic='blue';

ls{3}='-';
ls{2}='--';
ls{1}=':';
mi=1.25;
%ls{4}='-';

% %load data
cons=0;
for select=0:1
    AUC=zeros(10,total+1);
    figure('units','normalized','position',[0 0 1 1])
    for j=1:20
        if select==0
            filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
        else
            filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
        end
        load(filename)
        AUC(1,j+cons)=mean(powerP);
        AUC(2,j+cons)=mean(powerD);
        AUC(3,j+cons)=mean(powerM);
        AUC(4,j+cons)=mean(powerHHG);
        AUC(5,j+cons)=mean(powerMGC);
        AUC(6,j+cons)=mean(powerHSIC);
        AUC(7,j+cons)=mean(powerCorr);
        if select==0
        AUC(8,j+cons)=mean(powerSpearman);
        AUC(9,j+cons)=mean(powerKendall);
        AUC(10,j+cons)=mean(powerMIC);
        else
            AUC(8,j+cons)=mean(powerCCA);
        end
    end
    for j=1:total
        AUC(1:10,j)=AUC(1:10,j)./AUC(5,j);
    end
    AUC(1:10,total+1)=floor(mean(AUC(1:10,1:total-1),2)*100)/100;
    AUC=bsxfun(@min,AUC,mi-0.05);
    %
    x=1:total;
    hold on
    plot(x(1:total),AUC(5,1:total),'-.','MarkerSize',10,'LineWidth',3,'Color',MGC);
    for i=1:total
        text(x(i),AUC(1,i),'N','VerticalAlignment','middle','HorizontalAlignment','left','Color',mantel,'FontSize',fontSize-3);
        text(x(i),AUC(2,i),'D','VerticalAlignment','middle','HorizontalAlignment','left','Color',dcorr,'FontSize',fontSize-3);
        text(x(i),AUC(3,i),'C','VerticalAlignment','middle','HorizontalAlignment','left','Color',mcorr,'FontSize',fontSize-3);
        text(x(i),AUC(4,i),'G','VerticalAlignment','middle','HorizontalAlignment','left','Color',HHG,'FontSize',fontSize-3);
        text(x(i),AUC(6,i),'H','VerticalAlignment','middle','HorizontalAlignment','left','Color',hsic,'FontSize',fontSize-3);
        if select==0
        text(x(i),AUC(7,i),'P','VerticalAlignment','middle','HorizontalAlignment','left','Color',pcorr,'FontSize',fontSize-3);
        text(x(i),AUC(8,i),'S','VerticalAlignment','middle','HorizontalAlignment','left','Color',spearman,'FontSize',fontSize-3);
        text(x(i),AUC(9,i),'K','VerticalAlignment','middle','HorizontalAlignment','left','Color',kendall,'FontSize',fontSize-3);
        text(x(i),AUC(10,i),'I','VerticalAlignment','middle','HorizontalAlignment','left','Color',mic,'FontSize',fontSize-3);
        else
            text(x(i),AUC(7,i),'R','VerticalAlignment','middle','HorizontalAlignment','left','Color',pcorr,'FontSize',fontSize-3);
            text(x(i),AUC(8,i),'A','VerticalAlignment','middle','HorizontalAlignment','left','Color',spearman,'FontSize',fontSize-3);
        end
    end
    text(total+1,-0.07,'Average','VerticalAlignment','top','HorizontalAlignment','left','FontSize',fontSize);
    
    txt1=strcat('maNtel:',{' '},num2str(AUC(1,21)*100),'%');
    txt2=strcat('Dcorr:',{' '},num2str(AUC(2,21)*100),'%');
    txt3=strcat('mCorr:',{' '},num2str(AUC(3,21)*100),'%');
    txt4=strcat('hhG:',{' '},num2str(AUC(4,21)*100),'%');
    txt5=strcat('MGC:',{' '},num2str(AUC(5,21)*100),'%');
    txt6=strcat('Hsic:',{' '},num2str(AUC(6,21)*100),'%');
    if select==0
        txt7=strcat('Pearson:',{' '},num2str(AUC(7,21)*100),'%');
    txt8=strcat('Spearman:',{' '},num2str(AUC(8,21)*100),'%');
    txt9=strcat('Kendall:',{' '},num2str(AUC(9,21)*100),'%');
    txt10=strcat('mIc:',{' '},num2str(AUC(10,21)*100),'%');
    else
        txt7=strcat('RV:',{' '},num2str(AUC(7,21)*100),'%');
        txt8=strcat('ccA:',{' '},num2str(AUC(8,21)*100),'%');
    end
    adj=zeros(10,1);
     if select==0
        adj(1)=-0.025;
        adj(3)=0.05;
        adj(6)=-0.03;
        adj(8)=-0.05;
        %adj(7)=0.05;
        adj(9)=-0.1;
         adj(10)=-0.05;
         adj(4)=-0.02;
     else
         adj(1)=-0.08;
        adj(2)=-0.05;
        adj(8)=-0.04;
     end
    text(total+1,AUC(1,total+1)+adj(1),txt1,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mantel);
    text(total+1,AUC(2,total+1)+adj(2),txt2,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',dcorr);
    text(total+1,AUC(3,total+1)+adj(3),txt3,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mcorr);
    text(total+1,AUC(4,total+1)+adj(4),txt4,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',HHG);
    text(total+1,AUC(5,total+1)+adj(5),txt5,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',MGC);
    text(total+1,AUC(6,total+1)+adj(6),txt6,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',hsic);
    text(total+1,AUC(7,total+1)+adj(7),txt7,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',pcorr);
    text(total+1,AUC(8,total+1)+adj(8),txt8,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',spearman);
    if select==0
    text(total+1,AUC(9,total+1)+adj(9),txt9,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',kendall);
    text(total+1,AUC(10,total+1)+adj(10),txt10,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mic);
    end
%     h1=plot(x,AUC(1,:),'s','MarkerSize',10,'Color','red');
%     h2=plot(x,AUC(2,:),'o','MarkerSize',10,'Color','black');
%     h3=plot(x,AUC(3,:),'.','MarkerSize',20,'Color','magenta');
%     h4=plot(x,AUC(4,:),'*','MarkerSize',10,'Color',HHG);

    hold off
    xlim([0,total+1]);
    ylim([-0.05,mi]);
    set(gca,'FontSize',fontSize);
    set(gca,'XTick',[1,5,10,15,20],'FontSize',fontSize);
    set(gca,'YTick',0:0.25:mi,'YTickLabel',0:25:125,'FontSize',fontSize);
    
    xlabel('Simulation Type','FontSize',fontSize+5)%,...
    % 'Units', 'normalized','Position', [0.4, -0.18]);%, 'HorizontalAlignment', 'left')
    ylabel('Relative Power (%)','FontSize',fontSize+5)%, ...
    axis('square');
    if select~=1
        figNumber='1DPowerSummary';
        title('A. One-Dimensional Settings','FontSize',fontSize+7);
    else
        figNumber='HDPowerSummary';
        title('B. High-Dimensional Settings','FontSize',fontSize+7);
    end
    F.fname=strcat(pre2,figNumber);
    F.wh=[4.3 3]*2;
    print_fig(gcf,F)
end
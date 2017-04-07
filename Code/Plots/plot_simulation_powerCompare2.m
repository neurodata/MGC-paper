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
HHG   = [0.5,0.5,0.5];
mcorr='magenta';
mantel='red';
dcorr='blue';
MGC='cyan';

ls{3}='-';
ls{2}='--';
ls{1}=':';
%ls{4}='-';

% %load data
cons=0;
for select=0:1
    AUC=zeros(5,total+1);
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
        AUC(5,j+cons)=mean(powerMGCM);
    end
    for j=1:total
        AUC(1:5,j)=AUC(1:5,j)./AUC(5,j);
    end
    AUC(1:5,total+1)=floor(mean(AUC(1:5,1:total-1),2)*100)/100;
    
    %
    x=1:total;
    hold on
    plot(x(1:total),AUC(5,1:total),'-.','MarkerSize',10,'LineWidth',3,'Color',MGC);
    for i=1:total
        text(x(i),AUC(1,i),'A','VerticalAlignment','middle','HorizontalAlignment','left','Color',mantel,'FontSize',fontSize-3);
        text(x(i),AUC(2,i),'D','VerticalAlignment','middle','HorizontalAlignment','left','Color',dcorr,'FontSize',fontSize-3);
        text(x(i),AUC(3,i),'M','VerticalAlignment','middle','HorizontalAlignment','left','Color',mcorr,'FontSize',fontSize-3);
        text(x(i),AUC(4,i),'H','VerticalAlignment','middle','HorizontalAlignment','left','Color',HHG,'FontSize',fontSize-3);
    end
    text(total+1,AUC(5,total+1)+0.05,'Average','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',fontSize);
    
    txt1=strcat('mAntel:',{' '},num2str(AUC(1,21)*100),'%');
    txt2=strcat('Dcorr:',{' '},num2str(AUC(2,21)*100),'%');
    txt3=strcat('Mcorr:',{' '},num2str(AUC(3,21)*100),'%');
    txt4=strcat('Hhg:',{' '},num2str(AUC(4,21)*100),'%');
    txt5=strcat('MGC:',{' '},num2str(AUC(5,21)*100),'%');
    adj=zeros(5,1);
    if select==0
        adj(1)=-0.05;
        adj(2)=0.03;
    end
    text(total+1,AUC(1,total+1)+adj(1),txt1,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mantel);
    text(total+1,AUC(2,total+1)+adj(2),txt2,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',dcorr);
    text(total+1,AUC(3,total+1)+adj(3),txt3,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mcorr);
    text(total+1,AUC(4,total+1)+adj(4),txt4,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',HHG);
    text(total+1,AUC(5,total+1)+adj(5),txt5,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',MGC);
%     h1=plot(x,AUC(1,:),'s','MarkerSize',10,'Color','red');
%     h2=plot(x,AUC(2,:),'o','MarkerSize',10,'Color','black');
%     h3=plot(x,AUC(3,:),'.','MarkerSize',20,'Color','magenta');
%     h4=plot(x,AUC(4,:),'*','MarkerSize',10,'Color',HHG);

    hold off
    xlim([0,total+1]);
    ylim([-0.05,1.25]);
    set(gca,'FontSize',fontSize);
    set(gca,'XTick',[1,5,10,15,20],'FontSize',fontSize);
    set(gca,'YTick',0:0.25:1.25,'YTickLabel',0:25:125,'FontSize',fontSize);
    
    xlabel('Simulation Type','FontSize',fontSize+5)%,...
    % 'Units', 'normalized','Position', [0.4, -0.18]);%, 'HorizontalAlignment', 'left')
    ylabel('Relative Power (%)','FontSize',fontSize+5)%, ...
    axis('square');
    if select~=1
        figNumber='1DPowerSummary';
        title('One-Dimensional Settings','FontSize',fontSize+3);
    else
        figNumber='HDPowerSummary';
        title('High-Dimensional Settings','FontSize',fontSize+3);
    end
    F.fname=strcat(pre2,figNumber);
    F.wh=[3.8 3]*2;
    F.pdf=1;
    print_fig(gcf,F)
end
function []=plot_simulation_powerCompare3
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
thres=0.85;
mi=1000;
for select=0:1
    AUC=getData(select,total,pre1,thres);
    figure('units','normalized','position',[0 0 1 1])
    AUC(2,:)=AUC(3,:);
%     AUC(1:5,21)=floor(mean(AUC(1:5,1:19),2)*100)/100;
%     AUC=log(AUC);
    %
    x=1:total;
    hold on
    %plot(x(1:nn-1),AUC(5,1:nn-1),'-.','MarkerSize',10,'LineWidth',3,'Color',MGC);
    for i=1:total-1
        %text(x(i),AUC(1,i),'A','VerticalAlignment','middle','HorizontalAlignment','left','Color',mantel,'FontSize',fontSize-3);
        text(x(i),log10(min(AUC(2,i),mi)),'D','VerticalAlignment','middle','HorizontalAlignment','left','Color',dcorr,'FontSize',fontSize-3);
        %text(x(i),AUC(3,i),'M','VerticalAlignment','middle','HorizontalAlignment','left','Color',mcorr,'FontSize',fontSize-3);
        text(x(i),log10(min(AUC(4,i),mi)),'H','VerticalAlignment','middle','HorizontalAlignment','left','Color',HHG,'FontSize',fontSize-3);
        text(x(i),log10(min(AUC(5,i),mi)),'M','VerticalAlignment','middle','HorizontalAlignment','left','Color',MGC,'FontSize',fontSize-3);
    end
    text(total,log10(mi)+0.2,'Average','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',fontSize);
    
    txt1=strcat('mAntel:',{' '},num2str(floor(AUC(1,total))),'');
    txt2=strcat('Dcorr:',{' '},num2str(floor(AUC(2,total))),'');
    txt3=strcat('Mcorr:',{' '},num2str(floor(AUC(3,total))),'');
    txt4=strcat('Hhg:',{' '},num2str(floor(AUC(4,total))),'');
    txt5=strcat('MGC:',{' '},num2str(floor(AUC(5,total))),'');
    adj=zeros(5,1);
    if select==0
%         adj(1)=5;
%         adj(2)=0;
%         adj(3)=-5;
    end
    %text(nn,AUC(1,nn)+adj(1),txt1,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mantel);
    text(total,log10(min(AUC(2,total)+adj(2),mi)),txt2,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',dcorr);
    %text(nn,AUC(3,nn)+adj(3),txt3,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mcorr);
    text(total,log10(min(AUC(4,total)+adj(4),mi)),txt4,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',HHG);
    text(total,log10(min(AUC(5,total)+adj(5),mi)),txt5,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',MGC);
%     h1=plot(x,AUC(1,:),'s','MarkerSize',10,'Color','red');
%     h2=plot(x,AUC(2,:),'o','MarkerSize',10,'Color','black');
%     h3=plot(x,AUC(3,:),'.','MarkerSize',20,'Color','magenta');
%     h4=plot(x,AUC(4,:),'*','MarkerSize',10,'Color',HHG);

    hold off
    xlim([0,20]);
    ylim([0,3.5]);
    set(gca,'FontSize',fontSize);
    set(gca,'XTick',[1,5,10,15],'FontSize',fontSize);
    %ll=[{'0'};{'20'};{'40'};{'60'};{'80'};{'>=100'}];
    set(gca,'YTick',0:1:3,'YTickLabel',[1,10,100,1000], 'FontSize',fontSize);
    
    xlabel('Simulation Type','FontSize',fontSize+5)%,...
    % 'Units', 'normalized','Position', [0.4, -0.18]);%, 'HorizontalAlignment', 'left')
    ylabel(strcat('Sample Size to Achieve Power',{' '},num2str(thres)),'FontSize',fontSize+5)%, ...
    axis('square');
    if select~=1
        figNumber='1DPowerSummarySize';
        title('One-Dimensional Settings','FontSize',fontSize+3);
    else
        figNumber='HDPowerSummarySize';
        title('High-Dimensional Settings','FontSize',fontSize+3);
    end
    F.fname=strcat(pre2,figNumber);
    F.wh=[3.8 3]*2;
    F.pdf=1;
    print_fig(gcf,F)
end

function AUC=getData(select,total,pre1,thres)
AUC=zeros(5,total);
for j=1:total-1
    if select==0
        filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
        load(filename)
        tmp=find(powerP>thres,1);
        if isempty(tmp)
            tmp=100;
        end
        AUC(1,j)=tmp;
        
        tmp=find(powerD>thres,1);
        if isempty(tmp)
            tmp=100;
        end
        AUC(2,j)=tmp;
        
        tmp=find(powerM>thres,1);
        if isempty(tmp)
            tmp=100;
        end
        AUC(3,j)=tmp;
        
        tmp=find(powerHHG>thres,1);
        if isempty(tmp)
            tmp=100;
        end
        AUC(4,j)=tmp;
        
        tmp=find(powerMGCM>thres,1);
        if isempty(tmp)
            tmp=100;
        end
        AUC(5,j)=tmp;
        
        AUC(3,7)=120;
        AUC(3,8)=340;
        
        AUC(3,12)=400;
        AUC(3,13)=1000;
        
        AUC(3,14)=260;
        AUC(3,15)=500;
        AUC(3,16)=290;
        AUC(3,17)=300;
        AUC(3,18)=200;
        
        AUC(4,12)=330;
        AUC(4,13)=800;
        AUC(4,14)=120;
        AUC(4,15)=130;
        
        AUC(5,14)=120;
        AUC(5,15)=170;
    else
        AUC(3,1)=120;
        AUC(3,2)=340;
        AUC(3,3)=340;
        AUC(3,4)=340;
        AUC(3,5)=340;
        AUC(3,6)=340;
        AUC(3,7)=120;
        AUC(3,8)=340;
        AUC(3,9)=340;
        AUC(3,10)=340;
        AUC(3,11)=340;
        AUC(3,12)=340;
        AUC(3,13)=120;
        AUC(3,14)=340;
        AUC(3,15)=340;
        AUC(3,16)=340;
        AUC(3,17)=340;
        AUC(3,18)=340;
        AUC(3,19)=340;
        AUC(4,1)=120;
        AUC(4,2)=340;
        AUC(4,3)=340;
        AUC(4,4)=340;
        AUC(4,5)=340;
        AUC(4,6)=340;
        AUC(4,7)=120;
        AUC(4,8)=340;
        AUC(4,9)=340;
        AUC(4,10)=340;
        AUC(4,11)=340;
        AUC(4,12)=340;
        AUC(4,13)=120;
        AUC(4,14)=340;
        AUC(4,15)=340;
        AUC(4,16)=340;
        AUC(4,17)=340;
        AUC(4,18)=340;
        AUC(4,19)=340;
        AUC(5,1)=120;
        AUC(5,2)=340;
        AUC(5,3)=340;
        AUC(5,4)=340;
        AUC(5,5)=340;
        AUC(5,6)=340;
        AUC(5,7)=120;
        AUC(5,8)=340;
        AUC(5,9)=340;
        AUC(5,10)=340;
        AUC(5,11)=340;
        AUC(5,12)=340;
        AUC(5,13)=120;
        AUC(5,14)=340;
        AUC(5,15)=340;
        AUC(5,16)=340;
        AUC(5,17)=340;
        AUC(5,18)=340;
        AUC(5,19)=340;
    end
end

for j=1:size(AUC,1)
    AUC(j,total)=mean(AUC(j,1:total-1));%AUC(1:5,j)./AUC(5,j);
end
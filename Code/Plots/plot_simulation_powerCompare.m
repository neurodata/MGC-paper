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
%ls{4}='-';

% %load data
thres=0.85;
ind=1;
opt=1;
mi=500;
% if opt==1
%     mi=10;
% else
%     mi=1000;
% end
for select=0:1
    if select==0
        filename=strcat(pre1,'CorrIndSizeDim1.mat');
    else
        filename=strcat(pre1,'CorrIndSizeDimH.mat');
    end
    load(filename);
    AUC=SampleSize(2:end,:);
    AUC(:,total)=median(AUC(:,1:total-1),2);
%     if opt==1
%         AUC=AUC./repmat(AUC(ind,:),size(AUC,1),1);
%         AUC=floor(AUC*10)/10;
%     end
    figure('units','normalized','position',[0 0 1 1])
    x=1:total;
    hold on
%     AUC=AUC./AUC(ind,:);
%     AUC=floor(AUC*10)/10;
    for i=1:total-1       
        text(x(i),min(AUC(ind,i),mi),'M','VerticalAlignment','middle','HorizontalAlignment','left','Color',MGC,'FontSize',fontSize-3);     
        text(x(i),min(AUC(2,i),mi),'C','VerticalAlignment','middle','HorizontalAlignment','left','Color',mcorr,'FontSize',fontSize-3);      
        text(x(i),min(AUC(3,i),mi),'D','VerticalAlignment','middle','HorizontalAlignment','left','Color',dcorr,'FontSize',fontSize-3);
        text(x(i),min(AUC(4,i),mi),'N','VerticalAlignment','middle','HorizontalAlignment','left','Color',mantel,'FontSize',fontSize-3);
        text(x(i),min(AUC(5,i),mi),'G','VerticalAlignment','middle','HorizontalAlignment','left','Color',HHG,'FontSize',fontSize-3);
        text(x(i),min(AUC(6,i),mi),'H','VerticalAlignment','middle','HorizontalAlignment','left','Color',hsic,'FontSize',fontSize-3);
        if select==0
            text(x(i),min(AUC(7,i),mi),'P','VerticalAlignment','middle','HorizontalAlignment','left','Color',pcorr,'FontSize',fontSize-3);
            text(x(i),min(AUC(8,i),mi),'S','VerticalAlignment','middle','HorizontalAlignment','left','Color',spearman,'FontSize',fontSize-3);
            text(x(i),min(AUC(9,i),mi),'K','VerticalAlignment','middle','HorizontalAlignment','left','Color',kendall,'FontSize',fontSize-3);
            text(x(i),min(AUC(10,i),mi),'I','VerticalAlignment','middle','HorizontalAlignment','left','Color',mic,'FontSize',fontSize-3);
        else
            text(x(i),min(AUC(7,i),mi),'R','VerticalAlignment','middle','HorizontalAlignment','left','Color',pcorr,'FontSize',fontSize-3);
            text(x(i),min(AUC(8,i),mi),'A','VerticalAlignment','middle','HorizontalAlignment','left','Color',spearman,'FontSize',fontSize-3);
        end
    end
    text(total,-0.01*mi,'Median','VerticalAlignment','top','HorizontalAlignment','left','FontSize',fontSize);
    
    txt4=strcat('maNtel:',{' '},num2str(AUC(4,total)));
    txt3=strcat('Dcorr:',{' '},num2str(AUC(3,total)));
    txt2=strcat('mCorr:',{' '},num2str(AUC(2,total)));
    txt5=strcat('hhG:',{' '},num2str(AUC(5,total)));
    txt1=strcat('MGC:',{' '},num2str(AUC(ind,total)));
    txt6=strcat('Hsic:',{' '},num2str(AUC(6,total)));
    if select==0
        txt7=strcat('Pearson',{'>'},num2str(AUC(7,total)));
        txt8=strcat('Spearman',{'>'},num2str(AUC(8,total)));
        txt9=strcat('Kendall',{'>'},num2str(AUC(9,total)));
        txt10=strcat('mIc:',{' '},num2str(AUC(10,total)));
    else
        txt7=strcat('RV',{'>'},num2str(AUC(7,total)));
        txt8=strcat('ccA',{'>'},num2str(AUC(8,total)));
    end
    adj=zeros(10,1);
    if select==0
        adj(1)=-0.1*mi;
        adj(2)=0.1*mi;
        adj(3)=0.0*mi;
        adj(4)=0.05*mi;
        adj(5)=-0.05*mi;
        adj(6)=-0.025*mi;
        adj(7)=0*mi;
        adj(8)=-0.05*mi;
        adj(9)=-0.1*mi;
        adj(10)=0.075*mi;
    else
        adj(7)=0*mi;
        adj(8)=-0.05*mi;
%         adj(2)=-0.02*mi;
%         adj(5)=-0.03*mi;
%         adj(7)=-0.5*mi;
%         adj(8)=-0.55*mi;
    end
% %     %text(nn,AUC(1,nn)+adj(1),txt1,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mantel);
%     text(total,log10(min(AUC(3,total)+adj(3),mi)),txt3,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mcorr);
%     %text(nn,AUC(3,nn)+adj(3),txt3,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mcorr);
%     text(total,log10(min(AUC(4,total)+adj(4),mi)),txt4,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',HHG);
%     text(total,log10(min(AUC(5,total)+adj(5),mi)),txt5,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',MGC);
    text(total,min(AUC(4,total)+adj(4),mi),txt4,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mantel);
    text(total,min(AUC(3,total)+adj(3),mi),txt3,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',dcorr);
    text(total,min(AUC(2,total)+adj(2),mi),txt2,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mcorr);
    text(total,min(AUC(5,total)+adj(5),mi),txt5,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',HHG);
    text(total,min(AUC(ind,total)+adj(1),mi),txt1,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',MGC);
    text(total,min(AUC(6,total)+adj(6),mi),txt6,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',hsic);
    text(total,min(AUC(7,total)+adj(7),mi+adj(7)),txt7,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',pcorr);
    text(total,min(AUC(8,total)+adj(8),mi+adj(8)),txt8,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',spearman);
    if select==0
    text(total,min(AUC(9,total)+adj(9)+adj(9),mi+adj(9)),txt9,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',kendall);
    text(total,min(AUC(10,total)+adj(10),mi),txt10,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',fontSize,'Color',mic);
    end
%     h1=plot(x,AUC(1,:),'s','MarkerSize',10,'Color','red');
%     h2=plot(x,AUC(2,:),'o','MarkerSize',10,'Color','black');
%     h3=plot(x,AUC(3,:),'.','MarkerSize',20,'Color','magenta');
%     h4=plot(x,AUC(4,:),'*','MarkerSize',10,'Color',HHG);

    hold off
    xlim([0,20]);
     ylim([0,mi+10]);
%     ylim([0,10]);
    set(gca,'FontSize',fontSize);
    set(gca,'XTick',[1,5,10,15],'FontSize',fontSize);
    %ll=[{'0'};{'20'};{'40'};{'60'};{'80'};{'>=100'}];
    %set(gca,'YTick',0:1:3,'YTickLabel',[1,10,100,1000], 'FontSize',fontSize);
    
    xlabel('Simulation Type','FontSize',fontSize+5)%,...
    % 'Units', 'normalized','Position', [0.4, -0.18]);%, 'HorizontalAlignment', 'left')
    if select==0
    if opt==1
        ylabel(strcat('Relative Size to Achieve',{' '},num2str(thres*100),'% Power'),'FontSize',fontSize+5)%, ...
    else
        ylabel(strcat('Efficiency to Achieve',{' '},num2str(thres*100),'% Power'),'FontSize',fontSize+5)%, ...
    end
    end
    axis('square');
    if select~=1
        figNumber='1DPowerSummarySize';
        title('A. One-Dimensional Settings','FontSize',fontSize+7);
    else
        figNumber='HDPowerSummarySize';
        title('B. High-Dimensional Settings','FontSize',fontSize+7);
    end
    F.fname=strcat(pre2,figNumber);
    F.wh=[4.5 3]*2;
    print_fig(gcf,F)
end

% 
% % but can you also send me a normalized version?
% % that is, for each setting, divide by n for MGC?
% function AUC=getData(select,total,pre1,thres)
% AUC=zeros(5,total);
% for j=1:total-1
%     if select==0
%         filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
%         load(filename)
%         tmp=find(powerP>thres,1);
%         if isempty(tmp)
%             tmp=100;
%         end
%         AUC(1,j)=tmp;
%         
%         tmp=find(powerD>thres,1);
%         if isempty(tmp)
%             tmp=100;
%         end
%         AUC(2,j)=tmp;
%         
%         tmp=find(powerM>thres,1);
%         if isempty(tmp)
%             tmp=100;
%         end
%         AUC(3,j)=tmp;
%         
%         tmp=find(powerHHG>thres,1);
%         if isempty(tmp)
%             tmp=100;
%         end
%         AUC(4,j)=tmp;
%         
%         tmp=find(powerMGCM>thres,1);
%         if isempty(tmp)
%             tmp=100;
%         end
%         AUC(5,j)=tmp;
%         
%         AUC(3,7)=120;
%         AUC(3,8)=340;
%         
%         AUC(3,12)=400;
%         AUC(3,13)=1000;
%         
%         AUC(3,14)=260;
%         AUC(3,15)=500;
%         AUC(3,16)=290;
%         AUC(3,17)=300;
%         AUC(3,18)=200;
%         
%         AUC(4,12)=330;
%         AUC(4,13)=800;
%         AUC(4,14)=120;
%         AUC(4,15)=130;
%         
%         AUC(5,14)=120;
%         AUC(5,15)=170;
%     else
%         AUC(3,1)=20;
%         AUC(3,2)=20;
%         AUC(3,3)=25;
%         AUC(3,4)=180;
%         AUC(3,5)=130;      
%         AUC(3,6)=130;  
%         AUC(3,7)=170;   
%         AUC(3,8)=170;
%         AUC(3,9)=150;
%         AUC(3,10)=120;
%         AUC(3,11)=410;
%         AUC(3,12)=1000;
%         AUC(3,13)=1000;
%         AUC(3,14)=500;      
%         AUC(3,15)=680;
%         AUC(3,16)=1000;
%         AUC(3,17)=1000;
%         AUC(3,18)=400;
%         AUC(3,19)=500;
%         
%         AUC(4,1)=30;
%         AUC(4,2)=30;
%         AUC(4,3)=40;
%         AUC(4,4)=550;
%         AUC(4,5)=230;
%         AUC(4,6)=190;
%         AUC(4,7)=250;
%         AUC(4,8)=600;
%         AUC(4,9)=90;
%         AUC(4,10)=90;
%         AUC(4,11)=210;
%         AUC(4,12)=1000;
%         AUC(4,13)=1000;
%         AUC(4,14)=180;
%         AUC(4,15)=190;
%         AUC(4,16)=140;
%         AUC(4,17)=140;
%         AUC(4,18)=90;
%         AUC(4,19)=40;
%         
%         AUC(5,1)=20;
%         AUC(5,2)=20;
%         AUC(5,3)=25;
%         AUC(5,4)=180;
%         AUC(5,5)=130; 
%         AUC(5,6)=90;
%         AUC(5,7)=110;
%         AUC(5,8)=170;
%         AUC(5,9)=90;
%         AUC(5,10)=40;
%         AUC(5,11)=90;
%         AUC(5,12)=200;     
%         AUC(5,13)=1000;  
%         AUC(5,14)=140;        
%         AUC(5,15)=130;    
%         AUC(5,16)=100;
%         AUC(5,17)=100;
%         AUC(5,18)=80;
%         AUC(5,19)=130;
%     end
% end
% 
% for j=1:size(AUC,1)
%     AUC(j,total)=mean(AUC(j,1:total-1));%AUC(1:5,j)./AUC(5,j);
% end
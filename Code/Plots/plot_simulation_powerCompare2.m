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
fontSize=20;
%% Set colors
loca=[0,1,0];
glob= [0.5,0.5,0.5];
HHG   = [0.5,0.5,0.5];

ls{3}='-';
ls{2}='--';
ls{1}=':';
%ls{4}='-';

AUC=zeros(5,total);
% %load data
cons=0;
for select=0:1
    figure('units','normalized','position',[0 0 1 1])
    for j=1:total
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
    
    %
    x=1:total;
    hold on
    h1=plot(x,AUC(1,:),'+','MarkerSize',10,'Color','red');
    h2=plot(x,AUC(2,:),'o','MarkerSize',10,'Color','black');
    h3=plot(x,AUC(3,:),'.','MarkerSize',20,'Color','magenta');
    h4=plot(x,AUC(4,:),'x','MarkerSize',10,'Color',HHG);
    h5=plot(x,AUC(5,:),'-.','MarkerSize',10,'LineWidth',3,'Color','cyan');
    hold off
%         hold on
%         pv1=AUC(j,:);
%         %pv1=pv1(ss);
%         ord1=0.05*(1:total);
%         %     plot(pv2,ord1,'*','MarkerSize',8,'Color',gr);
%         
%         pv2=AUC(j+4,:);
%         %pv2=pv2(ss);
%         plot(pv2,ord1,'.','MarkerSize',20,'Color',bl);
        %plot(pv2,ord1,'*','MarkerSize',8,'Color',bl);
        %     if j==1
        %         h=legend('Oracle','Sample');
        %         set(h,'FontSize',fontSize, 'Units', 'normalized','Position',[0.06,0.75,0,0]);
        %         legend boxoff
        %     end
        %     plot(pv1,ord1,'.','MarkerSize',15,'Color',gr);
%         
%         plot(zeros(21,1),0:0.05:1,'--','Color',gr,'LineWidth',2)
%         mm1=mean(pv1);
%         mm2=mean(pv2);
%         pos=0;
%         %     plot(mm1,pos,'.','MarkerSize',30,'Color',gr);
%         %plot(mm2,pos,'*','MarkerSize',13,'Color',bl);
%         plot(mm2,pos,'.','MarkerSize',30,'Color',bl);
        
        %     a=text(mm1,pos+0.02,num2str(round(mm1*100)/100),'VerticalAlignment','bottom','HorizontalAlignment','left','Color',gr,'Interpreter','latex');
        %     set(a,'FontSize',fontSize);
        %     b=text(mm2,pos-0.03,num2str(round(mm2*100)/100),'VerticalAlignment','top','HorizontalAlignment','left','Color',bl,'Interpreter','latex');
        %     set(b,'FontSize',fontSize);
        
%         axis('square');
        hold off
        xlim([1,19]);
        ylim([0,1.25]);
        set(gca,'FontSize',fontSize);
%         set(gca,'XTick',[]); % Remove y axis ticks
        %set(gca,'XTickLabel',[0],'XTick',[0],'FontSize',fontSize);
%         labels={-1,[],[],[],0,[],[],[],1};
        set(gca,'XTick',[1,5,10,15,19],'FontSize',fontSize);
        set(gca,'YTick',0:0.25:1.25,'FontSize',fontSize);
        
        xlabel('Simulation Type','FontSize',fontSize)%,...
           % 'Units', 'normalized','Position', [0.4, -0.18]);%, 'HorizontalAlignment', 'left')
        ylabel('Mean Power Relative to Sample MGC','FontSize',fontSize)%, ...
        txt1=strcat('Mantel:',{' '},num2str(floor(mean(AUC(1,:))*100)/100));
        txt2=strcat('Dcorr:',{' '},num2str(floor(mean(AUC(2,:))*100)/100));
        txt3=strcat('Mcorr:',{' '},num2str(floor(mean(AUC(3,:))*100)/100));
        txt4=strcat('HHG:',{' '},num2str(floor(mean(AUC(4,:))*100)/100));
        txt5=strcat('MGC:',{' '},num2str(floor(mean(AUC(5,:))*100)/100));
        h=legend([h5,h1 h2 h3 h4],string(txt5{1}),string(txt1{1}),string(txt2{1}),string(txt3{1}),string(txt4{1}),'Location','NorthEastOutside');
           % 'Units', 'normalized','Position', [-0.18 0.38])%, 'HorizontalAlignment', 'left')
        
        %     switch j
        %         case 1
        %             tit=strcat('A. MGC - ',{' '}, str(j,:));
        %         case 2
        %             tit=strcat('B. MGC - ',{' '}, str(j,:));
        %         case 3
        %             tit=strcat('C. MGC - ',{' '}, str(j,:));
        %         case 4
        %             tit=strcat('D. MGC - ',{' '}, str(j,:));
        
%         txt1=strcat(tit,' Power(MGC)');
%         txt2=strcat({'     '},'- Power(',str(j,:),')');
%         txt2=txt2{1};
%         title({txt1,txt2},'FontSize',fontSize, ...
%             'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
%         hold off
    
    if select~=1
        figNumber='1DPowerSummary';
        title('One-Dimensional Settings');
    else
        figNumber='HDPowerSummary';
        title('High-Dimensional Settings');
    end
    F.fname=strcat(pre2,figNumber);
    F.wh=[4.5 3]*2;
    print_fig(gcf,F)
end
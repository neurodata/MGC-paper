function []=plot_simulation_permutation

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
pre2=strcat(rootDir,'Figures/Fig');% The folder to save figures

%% Set colors
%map2 = brewermap(128,'PRGn'); % brewmap
loca=[0,1,0];
glob= [1,0,1];
HHG   = [0.5,0.5,0.5];
mkSize=25;

ls{1}='-';
ls{2}='-';
ls{3}='--';

% set(groot,'defaultAxesColorOrder',map1);

%% 1-d
for i=1:2
    if i==1
        figNumber='1DPerm';
        filename=strcat(pre1,'CorrSimPermScale1-20Dim1');
    else
        figNumber='HDPerm';
        filename=strcat(pre1,'CorrSimPermScale1-20Dim2');
    end
    load(filename)
    
    x=1:20;
    a=2;
    if a==1
        ind=[1,2,3,7];
    else
        ind=[4,5,6,7];
    end
    p1=powerP(ind,:);
    p1=p1-repmat(p1(1,:),4,1);
    pp2=p1(2,:);
    pp1=p1(3,:);
    pp3=p1(4,:);
    [f1,xi1]=ksdensity(pp1,'support',[-1,1]);
    [f2,xi2]=ksdensity(pp2,'support',[-1,1]);
    [f3,xi3]=ksdensity(pp3,'support',[-1,1]);
    
    figure
    hold on
    h3=plot(xi3,f3,ls{3},'Color',HHG,'LineWidth',4);
    h2=plot(xi2,f2,ls{2},'Color',glob,'LineWidth',4);
    h1=plot(xi1,f1,ls{1},'Color',loca,'LineWidth',4);
    h=legend([h1 h2 h3],'True MGC','Mcorr','HHG','Location','NorthWest');
    plot(pp3,-0.1*ones(size(pp3)),'.','MarkerSize',mkSize,'Color',HHG);
    plot(pp2,-0.25*ones(size(pp2)),'.','MarkerSize',mkSize,'Color',glob);
    plot(pp1,-0.4*ones(size(pp1)),'.','MarkerSize',mkSize,'Color',loca);
    hold off
    set(gca,'FontSize',24);
    ylim([-0.5 4]);
    set(gca,'YTick',[]); % Remove y axis ticks
    
    legend boxoff
    % xlabel('Function Type','FontSize',24);
    %ylabel('Power Difference','FontSize',24);
    
    % xlim([-1,20]);
    % ylim([-0.6,1]);
    % set(gca,'XTick',[1,10,20]); % Remove x axis ticks
    %     set(gca,'YTick',[-0.5,0,0.5,1]); % Remove x axis ticks
    if i==1
        title('1-Dimensional Simulations','FontSize',24);
    else
        title('High-Dimensional Simulations','FontSize',24);
    end
    % grid on
    F.fname=strcat(pre2, figNumber);
    F.wh=[3 2.5]*2;
    print_fig(gcf,F)
    
    figure
hold on
p2=powerP(ind,:);
% plot(-1+ones(20,1)+randn(20,1)/10,p2(3,:),ls{5},'Color','k','LineWidth',4,'MarkerSize',mk);
plot(ones(20,1)+randn(20,1)/10,p2(1,:),ls{5},'Color',loca,'LineWidth',4,'MarkerSize',mk);
plot(1+ones(20,1)+randn(20,1)/10,p2(2,:),ls{5},'Color',glob,'LineWidth',4,'MarkerSize',mk);
plot(2+ones(20,1)+randn(20,1)/10,p2(4,:),ls{5},'Color',HHG,'LineWidth',4,'MarkerSize',mk);
hold off
% legend('True MGC','Estimated MGC','Mcorr','HHG','Location','NorthWest');
set(gca,'FontSize',24);
legend boxoff
% xlabel('Function Type','FontSize',24);
ylabel('Power','FontSize',24);
% xlim([1,20]);
% ylim([-0.6,1]);
set(gca,'XTick',[0,1,2,3],'XtickLabel',[{'True MGC'}; {'Estimated MGC'}; {'Mcorr'}; {'HHG'}]); % Remove x axis ticks
set(gca,'YTick',[0,0.5,1]); % Remove x axis ticks
ax=gca;
ax.XTickLabelRotation=45;
if i==1
        title('1-Dimensional Simulations','FontSize',24);
    else
        title('High-Dimensional Simulations','FontSize',24);
    end
grid on
xlim([-0.3, 3.3])
F.fname=strcat(pre2, figNumber, 'A');
F.wh=[3 2.5]*2;
print_fig(gcf,F)
end
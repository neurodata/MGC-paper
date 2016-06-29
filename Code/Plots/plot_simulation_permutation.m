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
mk=30;

ls{1}='.-';
ls{2}='.-';
ls{3}='.-';
ls{4}='.-';
ls{5}='.';
ls{6}='-';
 ind=[4,5,6,7];

% set(groot,'defaultAxesColorOrder',map1);

%% 1-d
for i=1:2
    if i==1
        figNumber='1DPerm';
        filename=strcat(pre1,'CorrSimPermScale1-20Dim1');
        tstr='1-Dimensional Simulations';
    else
        figNumber='HDPerm';
        filename=strcat(pre1,'CorrSimPermScale1-20Dim2');
        tstr='High-Dimensional Simulations';
    end
    load(filename)
    
    figure
    hold on
    p2=powerP(ind,:);
    p2=p2-repmat(p2(1,:),4,1);
    plot(-1+ones(20,1)+randn(20,1)/10,p2(3,:),ls{5},'Color',loca,'LineWidth',4,'MarkerSize',mk);
    %plot(ones(20,1)+randn(20,1)/10,p2(1,:),ls{5},'Color',loca,'LineWidth',4,'MarkerSize',mk);
    plot(0+ones(20,1)+randn(20,1)/10,p2(2,:),ls{5},'Color',glob,'LineWidth',4,'MarkerSize',mk);
    plot(1+ones(20,1)+randn(20,1)/10,p2(4,:),ls{5},'Color',HHG,'LineWidth',4,'MarkerSize',mk);
    hold off
    % legend('True MGC','Estimated MGC','Mcorr','HHG','Location','NorthWest');
    set(gca,'FontSize',24);
    legend boxoff
    % xlabel('Function Type','FontSize',24);
    ylabel('Power Difference','FontSize',24);
    % xlim([1,20]);
    % ylim([-0.6,1]);
    % set(gca,'XTick',[0,1,2,3],'XtickLabel',[{'True MGC'}; {'Estimated MGC'}; {'Mcorr'}; {'HHG'}]); % Remove x axis ticks
    set(gca,'XTick',[0,1,2],'XtickLabel',[{'True MGC'}; {'Mcorr'}; {'HHG'}]); % Remove x axis ticks
    set(gca,'YTick',[-1,-0.5,0,0.5,1]); % Remove x axis ticks
    ax=gca;
    ax.XTickLabelRotation=45;
    title(tstr,'FontSize',24);
    grid on
    xlim([-0.3, 2.3])
    F.fname=strcat(pre2, figNumber);
    F.wh=[3 2.5]*2;
    print_fig(gcf,F)
end
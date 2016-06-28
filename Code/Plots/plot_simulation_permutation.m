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
mkSize=20;

ls{1}='-';
ls{2}='-';
ls{3}='--';

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
    legend([h1 h2 h3],'True MGC','Mcorr','HHG','Location','NorthWest');
    
    for j=1:3
        switch j
            case 1
                pv=pp3;
                cc=HHG;
                c=-0.1;
            case 2
                pv=pp2;
                cc=glob;
                c=-0.3;
            case 3
                pv=pp1;
                cc=loca;
                c=-0.5;
        end
        pv=sort(pv,'ascend');
        ord=0*ones(size(pv));
%         for jj=2:length(pv);
%             if pv(jj)-pv(jj-1)<0.001
%                 ord(jj)=ord(jj-1)+0.1;
%                 %pv(jj)=pv(jj-1);
%             end
%         end
        
        plot(pv,c*ones(size(pv))-ord,'.','MarkerSize',mkSize,'Color',cc);
        
    end
    hold off
    set(gca,'FontSize',24);
    ylim([-0.6 4]);
    set(gca,'YTick',[]); % Remove y axis ticks
    
    legend boxoff
    % xlabel('Function Type','FontSize',24);
    %ylabel('Power Difference','FontSize',24);
    
    % xlim([-1,20]);
    % ylim([-0.6,1]);
    % set(gca,'XTick',[1,10,20]); % Remove x axis ticks
    %     set(gca,'YTick',[-0.5,0,0.5,1]); % Remove x axis ticks
        title(tstr,'FontSize',24);
    % grid on
    F.fname=strcat(pre2, figNumber);
    F.wh=[3 2.5]*2;
    print_fig(gcf,F)
end
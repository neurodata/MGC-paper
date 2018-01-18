function []=plot_simulation_heatmaps2
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
repP=1000;
map2 = brewermap(128,'BuPu'); % brewmap
glob= [0.5,0.5,0.5];

figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
lw=1;
figNumber='1DHeat';
nn=60;
rep=100;
load(strcat(pre1,'CorrIndTest1DHeat.mat'));
thres=0.025;
for j=1:total
    subplot(s,t,j)
    nn=60;
    ph=zeros(nn,nn);
    for r=1:rep
        [x,y]=CorrSampleGenerator(j,nn,1,1,1);
        [~,tmp,~]=MGCSampleStat(x,y);
        ph=ph+tmp/rep;
    end
    ph=ph';
    titlechar=strcat(num2str(j),'.',{' '},CorrSimuTitle(j));
    kmin=2;
    %     thres=0.8;
    % ind=[find(max(power2,[],1)>=thres,1) lim];
    % lim=min(ind);
    %     ind=find(numRange==nn);
    %     if isempty(ind)
    %         ind=1;
    %     end
    %     ph=powerMLocal(kmin:numRange(ind),kmin:numRange(ind),ind)';
    %     tt=find(sum(ph,2)==0,1,'first');
    %     if isempty(tt)==false && tt~=1;
    %         ph(tt:end,:)=repmat(ph(tt-1,:),numRange(ind)-tt,1);
    %     end
    %     tt=find(sum(ph,1)==0,1,'first');
    %     if isempty(tt)==false && tt~=1;
    %         ph(:,tt:end)=repmat(ph(:,tt-1),1,numRange(ind)-tt);
    %     end
    hold on
    imagesc(ph);
    set(gca,'YDir','normal')
    colormap(map2)
    tmax=max(max(max(ph)),thres);
    caxis([tmax/2 tmax]);
    
    optimalInd=find(ph==tmax,1,'last');
    if ~isempty(optimalInd)
        [J,I]=ind2sub(size(ph),optimalInd);
        %         Ymin=min(I)-1;
        %         Ymax=max(I)-1;
        %         Xmin=min(J)-1;
        %         Xmax=max(J)-1;
        
        if j<total
            %             plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
            %             plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
            %             plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
            %             plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
            plot(I-1,J-1,'go','markerSize',6,'linewidth',3)
        end
    end
    %plot(nn-1,nn-1,'.','markerSize',24,'MarkerFaceColor',glob,'Color',glob)
    hold off
    xlim([1,nn-1]);
    ylim([1,nn-1]);
    set(gca,'FontSize',14);
    set(gca,'XTick',[1,round(nn/2)-1,nn-1],'XTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
    set(gca,'YTick',[1,round(nn/2)-1,nn-1],'YTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
    title(titlechar,'FontSize',14);
    if j~=1
        set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
        set(gca,'XColor','none')
        set(gca,'YColor','none')
    end
    axis('square')
end
xlabel('# X Neighbors','position',[-170 -12],'FontSize',24);
ylabel('# Y Neighbors','position',[-423 165],'FontSize',24);

%colorbar
% h=colorbar('Ticks',[0,thres/2,thres]);
% pos2=get(gca,'pos');
% pos=get(h,'pos');
% set(h,'position',[pos(1)+0.05 pos(2)-0.015 pos(3)+0.003 pos(4)+0.02]);
% set(gca,'position',[pos2(1) pos2(2) pos2(3) pos2(4)]);
% set(h,'FontSize',14);
h=suptitle(strcat('MGC Image'));
set(h,'FontSize',26,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)

figure('units','normalized','position',[0 0 1 1])
figNumber='HDHeat';
load(strcat(pre1,'CorrIndTestHDHeat.mat'));
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    thres2=0.5;
    ind=find(max(powerMGCM,[],1)>=thres2,1,'last');
    if isempty(ind)
        ind=1;
    end
    subplot(s,t,j)
    nn=100;
    ph=zeros(nn,nn);
    for r=1:rep
        [x,y]=CorrSampleGenerator(j,nn,dimRange(ind),1,0);
        [~,tmp,~]=MGCSampleStat(x,y);
        ph=ph+tmp/rep;
    end
    ph=ph';
    titlechar=strcat(num2str(j),'.',{' '},CorrSimuTitle(j));
    kmin=2;
    %     thres=0.8;
    % ind=[find(max(power2,[],1)>=thres,1) lim];
    % lim=min(ind);
    %     ind=find(numRange==nn);
    %     if isempty(ind)
    %         ind=1;
    %     end
    %     ph=powerMLocal(kmin:numRange(ind),kmin:numRange(ind),ind)';
    %     tt=find(sum(ph,2)==0,1,'first');
    %     if isempty(tt)==false && tt~=1;
    %         ph(tt:end,:)=repmat(ph(tt-1,:),numRange(ind)-tt,1);
    %     end
    %     tt=find(sum(ph,1)==0,1,'first');
    %     if isempty(tt)==false && tt~=1;
    %         ph(:,tt:end)=repmat(ph(:,tt-1),1,numRange(ind)-tt);
    %     end
    hold on
    imagesc(ph);
    set(gca,'YDir','normal')
    colormap(map2)
    tmax=max(max(max(ph)),thres);
    caxis([tmax/2 tmax]);
    
    optimalInd=find(ph==tmax);
    if ~isempty(optimalInd)
        if find(optimalInd==nn*nn)
            optimalInd=nn*nn;
        else
            optimalInd=optimalInd(end);
        end
        [J,I]=ind2sub(size(ph),optimalInd);
        %         Ymin=min(I)-1;
        %         Ymax=max(I)-1;
        %         Xmin=min(J)-1;
        %         Xmax=max(J)-1;
        
        if j<total
            %             plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
            %             plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
            %             plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
            %             plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
            plot(I-1,J-1,'go','markerSize',6,'linewidth',3)
        end
    end
    hold off
    xlim([1,nn-1]);
    ylim([1,nn-1]);
    set(gca,'FontSize',14);
    set(gca,'XTick',[1,round(nn/2)-1,nn-1],'XTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
    set(gca,'YTick',[1,round(nn/2)-1,nn-1],'YTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
    title(titlechar,'FontSize',14);
    if j~=1
        set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
        set(gca,'XColor','none')
        set(gca,'YColor','none')
    end
    axis('square')
end
xlabel('# X Neighbors','position',[-290 -20],'FontSize',24);
ylabel('# Y Neighbors','position',[-715 277],'FontSize',24);

% h=colorbar('Ticks',[0,thres/2,thres]);
% pos2=get(gca,'pos');
% pos=get(h,'pos');
% set(h,'position',[pos(1)+0.05 pos(2)-0.015 pos(3)+0.003 pos(4)+0.02]);
% set(gca,'position',[pos2(1) pos2(2) pos2(3) pos2(4)]);
% axis('square');

% set(h,'FontSize',14);
h=suptitle(strcat('MGC Image'));
set(h,'FontSize',26,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)
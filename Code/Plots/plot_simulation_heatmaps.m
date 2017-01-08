function []=plot_simulation_heatmaps(figNumber)
% Used to plot figure 1-8 used in tex. Run like

% total is usually 20.
% pre1 specifies the location to load data.
% pre2 specifies the location to save pictures.

%%% File path searching
if nargin<1
    figNumber='HDHeat';
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
pre2=strcat(rootDir,'Figures/Fig');% The folder to save figures

total=20;
repP=100;
map2 = brewermap(128,'BuPu'); % brewmap

figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
lw=1;
if strcmp(figNumber,'1DHeat')
    nn=60;
    
    for j=1:total
        filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
        load(filename)
        subplot(s,t,j)
        titlechar=strcat(num2str(j),'.',{' '},CorrSimuTitle(j));
        kmin=2;
        thres=0.8;
        % ind=[find(max(power2,[],1)>=thres,1) lim];
        % lim=min(ind);
        ind=find(numRange==nn);
        if isempty(ind)
            ind=1;
        end
        ph=powerMLocal(kmin:numRange(ind),kmin:numRange(ind),ind)';
        tt=find(sum(ph,2)==0,1,'first');
        if isempty(tt)==false && tt~=1;
            ph(tt:end,:)=repmat(ph(tt-1,:),numRange(ind)-tt,1);
        end
        tt=find(sum(ph,1)==0,1,'first');
        if isempty(tt)==false && tt~=1;
            ph(:,tt:end)=repmat(ph(:,tt-1),1,numRange(ind)-tt);
        end
        hold on
        imagesc(ph);
        set(gca,'YDir','normal')
        colormap(map2)
        caxis([0 thres])
        
        if j<total
            %             phmax=max(max(ph));
            testRepeat=1;
            while testRepeat==1;
                [x,y]=CorrSampleGenerator(j,nn,1,1,1);
                [pMGC,statMGC,pLocalCorr,localCorr,optimalInd]=MGCPermutationTest(squareform(pdist(x)),squareform(pdist(y)),repP,'mcor');
                if pMGC<0.05 && ((sum(optimalInd==n^2)==0 && j>5) || (sum(optimalInd==n^2)==1 && j<=5))
                    testRepeat=0;
                end
            end
            
            %     [~,~,~,optimalInd]=FindLargestRectangles((ph>=phmax-powerThres), [0 0 1]);
            %     optimalInd=find(optimalInd==1);
            [J,I]=ind2sub(size(localCorr),optimalInd);
            Ymin=min(I)-1;
            Ymax=max(I)-1;
            Xmin=min(J)-1;
            Xmax=max(J)-1;
            
            plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
            plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
            plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
            plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
            localCorr(Xmin:Xmax,Ymin:Ymax)=localCorr(Xmin:Xmax,Ymin:Ymax)+10;
            statMGC=statMGC+10;
            op2=find(localCorr==statMGC);
            if isempty(op2)
                k=Xmax+1;
                l=Ymax+1;
            else
                [k,l]=ind2sub(size(localCorr),op2);
            end
            plot(k-1,l-1,'go','markerSize',5)
        end
        hold off
        xlim([1,nn-1]);
        ylim([1,nn-1]);
        set(gca,'FontSize',14);
        set(gca,'XTick',[1,round(nn/2)-1,nn-1],'XTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
        set(gca,'YTick',[1,round(nn/2)-1,nn-1],'YTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
        title(titlechar,'FontSize',14, ...
            'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
        if j~=1
            set(gca,'XTick',[]); % Remove x axis ticks
            set(gca,'YTick',[]); % Remove y axis ticks
        end
        axis('square')
    end
    xlabel('# X Neighbors','position',[-200 -12],'FontSize',24);
    ylabel('# Y Neighbors','position',[-500 180],'FontSize',24);
    
    %colorbar
    h=colorbar('Ticks',[0,thres/2,thres]);
    set(h,'FontSize',14);
    h=suptitle(strcat('One-Dimensional Multiscale Power Maps'));
else
    figNumber='HDHeat';
    for j=1:total
        filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
        load(filename)
        subplot(s,t,j)
        titlechar=strcat(num2str(j),'.',{' '},CorrSimuTitle(j));
        kmin=2;thres=0.5;
        ind=find(max(powerMGCM,[],1)>=thres,1,'last');
        if isempty(ind)
            ind=1;
        end
        ph=powerMLocal(kmin:n,kmin:n,ind)';
        tt=find(sum(ph,2)==0,1,'first');
        if isempty(tt)==false && tt~=1;
            ph(tt:end,:)=repmat(ph(tt-1,:),n-tt,1);
        end
        tt=find(sum(ph,1)==0,1,'first');
        if isempty(tt)==false && tt~=1;
            ph(:,tt:end)=repmat(ph(:,tt-1),1,n-tt);
        end
        hold on
        imagesc(ph);
        set(gca,'YDir','normal')
        colormap(map2)
        caxis([0 thres])
        
        if j<total
            %phmax=max(max(ph));
            testRepeat=1;
            while testRepeat==1;
                [x,y]=CorrSampleGenerator(j,n,dimRange(ind),1,0);
                [pMGC,statMGC,pLocalCorr,localCorr,optimalInd]=MGCPermutationTest(squareform(pdist(x)),squareform(pdist(y)),repP,'mcor');
                if pMGC<0.05 && ((sum(optimalInd==n^2)==0 && j>5) || (sum(optimalInd==n^2)==1 && j<=5))
                    testRepeat=0;
                end
            end
            
            %     [~,~,~,optimalInd]=FindLargestRectangles((ph>=phmax-powerThres), [0 0 1]);
            %     optimalInd=find(optimalInd==1);
            [J,I]=ind2sub(size(localCorr),optimalInd);
            Ymin=min(I)-1;
            Ymax=max(I)-1;
            Xmin=min(J)-1;
            Xmax=max(J)-1;
            
            plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
            plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
            plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
            plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
            localCorr(Xmin:Xmax,Ymin:Ymax)=localCorr(Xmin:Xmax,Ymin:Ymax)+10;
            statMGC=statMGC+10;
            op2=find(localCorr==statMGC);
            if isempty(op2)
                k=Xmax+1;
                l=Ymax+1;
            else
                [k,l]=ind2sub(size(localCorr),op2);
            end
            plot(k-1,l-1,'go','markerSize',5,'linewidth',8)
        end
        hold off
        xlim([1,n-1]);
        ylim([1,n-1]);
        set(gca,'XTick',[1,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n]); % Remove x axis ticks
        set(gca,'YTick',[1,round(n/2)-1,n-1],'YTickLabel',[2,round(n/2),n]); % Remove x axis ticks
        if j~=1
            set(gca,'XTick',[]); % Remove x axis ticks
            set(gca,'YTick',[]); % Remove y axis ticks
        end
        set(gca,'FontSize',14);
        title(titlechar,'FontSize',14);
        axis('square');
    end
    xlabel('# X Neighbors','position',[-340 -20],'FontSize',24);
    ylabel('# Y Neighbors','position',[-840 300],'FontSize',24);
    
    h=colorbar('Ticks',[0,thres/2,thres]);
    set(h,'FontSize',14);
    h=suptitle(strcat('High-Dimensional Multiscale Power Maps'));
end
set(h,'FontSize',24,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 5]*2;
print_fig(gcf,F)
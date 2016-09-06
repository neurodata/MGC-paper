function [pMGC]=plot_simulation_pValueAndPower(type,n,dim,noise)
% Compare p-value heatmap with power heatmap
if nargin<1
    type=13;
end
if nargin<2
    n=100;
end
if nargin<3
    dim=1;
end
if nargin<4
    noise=1;
end
rep=200;
fontSize=20;

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
map2 = brewermap(128,'GnBu'); % brewmap

% For p-value:
% Get distance matrix
[x, y]=CorrSampleGenerator(type,n,dim,1, noise);
C=squareform(pdist(x));
D=squareform(pdist(y));

% Calculate permutation p-value
tA=LocalCorr(C,D,2);
test=SampleMGC(tA);
tN=zeros(rep,n,n);
pMLocal=zeros(n,n);
pMGC=0;
for r=1:rep;
    per=randperm(n);
    tmp=LocalCorr(C,D(per,per),2);
    tmp2=SampleMGC(tmp);
    [t1,t2]=size(tmp);
    tN(r,1:t1,1:t2)=tmp;
    pMLocal=pMLocal+(tmp>=tA)/rep;
    pMGC=pMGC+(tmp2>=test)/rep;
end
if t1<n
    %tN(r,t1+1:n,1:t2)=repmat(tmp(t1,:),n-t1,1);
    pMLocal(t1+1:n,1:t2)=repmat(pMLocal(t1,:),n-t1,1);
end
if t2<n
    %tN(r,:,t2+1:n)=repmat(tmp(r,:,t2),1,n-t2);
    pMLocal(:,t2+1:n)=repmat(pMLocal(:,t2),n-t1,1);
end
if pMGC==0 || min(min(pMLocal(2:end,2:end)))==0
    pMLocal=pMLocal+1/rep;
    pMGC=pMGC+1/rep;
end
pMLocal(pMLocal>1)=1;

h=figure(type);
set(h,'units','normalized','position',[0 0 1 1]);
ax=subplot(2,2,1);
ph=pMLocal(2:end,2:end)';
% ph(ph<=eps)=0.0005;
imagesc(log(ph)); %log(ph)-min(log(ph(:))));
axis('square')
set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
cmap=map2;
colormap(ax,flipud(cmap));
% caxis([0 1]);
cticks=[0.001, 0.01, 0.1, 0.5];
h=colorbar('Ticks',log(cticks),'TickLabels',cticks);%,'location','westoutside');
set(h,'FontSize',fontSize);
set(gca,'XTick',[1,round(n/2)-1,n-1],'YTick',[1,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',16);
%xlabel('# of Neighbors for X','FontSize',16)
%ylabel('# of Neighbors for Y','FontSize',16) %,'Rotation',0,'position',[-7,20]);
xlim([1 n-1]);
ylim([1 n-1]);
title('Multiscale P-value Map')

%% For Power
[~,~,powerMLocal]=CorrIndTest(type,n,dim,1,rep, rep,noise,0.05,[0,1,0,0]);
kmin=2;
ph=powerMLocal(kmin:end,kmin:end)';
tt=find(sum(ph,2)==0,1,'first');
if isempty(tt)==false && tt~=1;
    ph(tt:end,:)=repmat(ph(tt-1,:),n-tt,1);
end
tt=find(sum(ph,1)==0,1,'first');
if isempty(tt)==false && tt~=1;
    ph(:,tt:end)=repmat(ph(:,tt-1),1,n-tt);
end
ax=subplot(2,2,2);
imagesc(ph);
axis('square')
set(gca,'YDir','normal')
colormap(map2)
caxis([0 1])
cticks=[0,0.5,1];
h=colorbar('Ticks',cticks);%,'location','westoutside');
set(gca,'FontSize',fontSize);
nn=size(ph,1)+1;
set(gca,'XTick',[1,round(nn/2)-1,nn-1],'XTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
set(gca,'YTick',[1,round(nn/2)-1,nn-1],'YTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
xlim([1 nn-1]);
ylim([1 nn-1]);
title('Multiscale Power Map')
h=suptitle(CorrSimuTitle(type));
set(h,'FontSize',26);% for 1-Dimensional Simulations'));
pos = get(ax,'position');

%Optimal Scale map
pMLocal=pMLocal';
% ind=find(pAll==min(min(pAll(2:end,2:end))));
if min(min(pMLocal))>pMGC
    pMGC=min(min(pMLocal));
end
% % [pval,phS]=MGCScaleVerify(pMLocal);
% phS=(pMLocal<=pMGC);
[~,~,~,phS]=FindLargestRectangles((pMLocal<=pMGC), [0 0 1],[2,2]);
% phS=find(phS==1);
if pMLocal(end)<=pMGC || sum(sum(phS))==0
    phS(end)=1;
end
% ind
% phS=zeros(size(pAll));
% phS(ind)=1;
ax=subplot(2,2,3);
imagesc(phS);
axis('square')
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize);
nn=size(phS,1);
set(gca,'XTick',[1,round(nn/2)-1,nn-1],'XTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
set(gca,'YTick',[1,round(nn/2)-1,nn-1],'YTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
xlim([1 nn-1]);
ylim([1 nn-1]);
title('Optimal Scale by P-value')
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

% pval=unique(pAll(ind));


%Optimal Scale map
ind=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
phS=zeros(size(ph));
phS(ind)=1;
ax=subplot(2,2,4);
imagesc(phS);
axis('square')
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize);
nn=size(ph,1)+1;
set(gca,'XTick',[1,round(nn/2)-1,nn-1],'XTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
set(gca,'YTick',[1,round(nn/2)-1,nn-1],'YTickLabel',[2,round(nn/2),nn]); % Remove x axis ticks
xlim([1 nn-1]);
ylim([1 nn-1]);
title('Optimal Scale by Power')
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);


%%
F.fname=strcat(rootDir, 'Figures/Auxiliary/PowerEst_type',num2str(type),'_n', num2str(n),'_d', num2str(dim));
F.wh=[8 5]*2;
print_fig(gcf,F)
function []=plot_realData3(opt)
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
if nargin<1
    opt=1;
end
pre2=strcat(rootDir,'Figures/FigReal',num2str(opt));% The folder to save figures

%% figure stuff
fontSize=9;
s=3;
t=4;
cmap=zeros(2,3);
map2 = brewermap(128,'PiYG'); % brewmap
gr=map2(120,:);
pu=map2(8,:);
cmap(1,:) = pu;
cmap(2,:) = gr;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

fnames={'CorrPermDistTestTypeBrainCxP.mat'; ...
    'CorrPermDistTestTypeBrainLMRxY.mat'; ...
    'CorrPermDistTestTypeMigrainxCCI.mat'};
fnames2={'BrainCP.mat'; ...
    'BrainHippoShape.mat'; ...
    'Semipar.mat'};
% xlabs={ '# Activity Neighbors'; ...
%         '# Shape Neighbors';...
%         '# Graph Neighbors'};
% ylabs={ '# Personality Neighbors'; ...
%         '# Disease Neighbors';...
%         '# Creativity Neighbors'};
% tits= {'A. Brain Activity vs. Personality'; ...
%     'B. Brain Shape vs. Disorder';...
%     'C. Brain Graph vs. Creativity'};
% tits={'A';'B';''};


%% loop maps
figure(1);
filename=strcat(pre1,fnames{opt});
load(filename);
glob= [0.5,0.5,0.5];
    
load(strcat(rootDir,'Data/Preprocessed/',fnames2{opt}))
switch opt
    case 1
        C=distC;
        D=distP;   
        per=1:size(distC,1);
    case 2
        C=LMRS;
        D=squareform(pdist((Label+unifrnd(0,0.01,n,1))));
        [~,per]=sort(Label);
    case 3
        C=distMigrain(ind,ind);
        D=squareform(pdist(cci));
        D=D(ind,ind);
        [~,per]=sort(cci(ind));
end
C=C(per,per)/max(max(C));
D=D(per,per)/max(max(D));
n=size(C,1);
% A Mcorr
H=eye(n)-(1/n)*ones(n,n);
A=H*C-C/n;
B=D*H-D/n;
for j=1:n
    A(j,j)=0;
    B(j,j)=0;
end
mcorrH=A.*B;
sameBar=1;
if sameBar==1
    minC=floor(min(min([A,B]))*10)/10;maxC=ceil(max(max([A,B]))*10)/10;
else
    minC=floor(min(min([A]))*10)/10;maxC=ceil(max(max([A]))*10)/10;
end
minC=min(minC,-maxC);maxC=max(maxC,-minC);

neighbor=optimalInd;
k=ceil(neighbor/n);
l=neighbor-n*(k-1);
[~,tind]=max(k.*l);
k=k(tind);
l=l(tind);

RC=DistRanks(C);
RD=DistRanks(D)';
RC=(RC>l);
RD=(RD>k);
A_MGC=A;
B_MGC=B;
A_MGC(RC)=0;
B_MGC(RD)=0;
% A
C_MGC=(A_MGC).*(B_MGC);

MH=max(max(C_MGC(2:end,2:end)));
mH=min(min(C_MGC(2:end,2:end)));
mH=min(mH,-MH);MH=max(MH,-mH);

% A Mantel
% ax=subplot(s,t,2); %
% ax=subplot('Position',[left(2), bottom(3), width, height]);
% ax=subplot(3,3,2);
% i=opt;
ax=subplot(s,t,1);
n=size(C,1);
hold on
sameBar=0;
colormap(map2);
% if sameBar==1
%     maxC=ceil(max(max([C,D]))*10)/10;
%     minC=ceil(min(min([C,D]))*10)/10;
% else
%     maxC=ceil(max(max(C))*10)/10;
%     minC=ceil(min(min(C))*10)/10;
% end
% minC=min(minC,-maxC);maxC=max(maxC,-minC);
imagesc(C');
% colormap(cmap);
caxis([minC,maxC]);

set(gca,'YDir','normal')
set(gca,'FontSize',fontSize); % Remove y axis ticks
% xlabel('# X Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
%     'Position', [0 -0.16], 'HorizontalAlignment', 'left')
% ylabel('# Y Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
%     'Position', [-0.25 0.4], 'VerticalAlignment', 'bottom')
% 'Position', [-0.2 -0.05], 'HorizontalAlignment', 'Left')
% text(45,110,'$\tilde{A}$','interpreter','latex','FontSize',fontSize)
% title([{'1. Mantel'}; {'(pairwise distances)'}; {' ' }],'FontSize',fontSize-1, 'Units', 'normalized', ...
%     'Position', [0 1], 'HorizontalAlignment', 'left');
title('$\tilde{A}$','interpreter','latex','FontSize',fontSize);
% title([{'Mantel'}; {' '}],'FontSize',fontSize);
%clean_panel(ax,map2,pos,id,n,col,fontSize)
set(gca,'visible','on')
axis('square');
% set(gca,'XTick',[2,round(n/2),n],'XtickLabel',[1,round(n/2),n]); % Remove x axis ticks
% set(gca,'YTick',[2.5,round(n/2),n],'YtickLabel',[1,round(n/2),n]); % Remove x axis ticks
xlim([1,n-1]);
ylim([1,n-1]);
hold off

% B Mantel
% ax=subplot(s,t,t+2);
ax=subplot(s,t,t+1);
hold all
% if sameBar~=1
%     maxC=ceil(max(max(D))*10)/10;
% end
imagesc(D');
set(gca,'FontSize',16)
% colormap(ax,map3);
caxis([minC,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize); % Remove y axis ticks
title('$\tilde{B}$','interpreter','latex','FontSize',fontSize);
axis('square');
%clean_panel(ax,map2,pos,id,n,col,fontSize)
xlim([1,n-1]);
ylim([1,n-1]);

% C Mantel
% ax=subplot(s,t,2*t+2);
ax=subplot(s,t,2*t+1);
hold all
n=size(C,1);
mantelH=(C-sum(sum(C))/n/(n-1)).*(D-sum(sum(D))/n/(n-1));
% MH=ceil(max(max(mantelH(2:end,2:end))));
% mH=floor(min(min(mantelH(2:end,2:end))));
% mH=min(mH,-MH);
% MH=max(MH,-mH);
imagesc(mantelH');
xlim([1,n-1]);
ylim([1,n-1]);
set(gca,'YDir','normal')
title('$$\tilde{C} = \tilde{A} \circ \tilde{B}$$','FontSize',fontSize,'interpreter','latex');
caxis([mH,MH]);
axis('square');
% colorbar('location','westoutside')
%clean_panel(ax,map2,pos,id,n,col,fontSize)


%% Mcorr

ax=subplot(s,t,2);
%ax=subplot('Position',[left(3), bottom(3), width, height]);
hold all
imagesc(A');
set(gca,'YDir','normal')
caxis([minC,maxC]);
title('$A$','interpreter','latex','FontSize',fontSize)
% title([{'2. Mcorr'}; {'(single center)'}; {' '}],'FontSize',fontSize-1, 'Units', 'normalized', ...
%     'Position', [0 1], 'HorizontalAlignment', 'left')
axis('square');
xlim([1,n-1]);
ylim([1,n-1]);
% clean_panel(ax,map2,pos,id,n,col,fontSize)


% B MCorr
% ax=subplot(s,t,t+3);
ax=subplot(s,t,t+2);
hold all
% if sameBar==1
%     minD=floor(min(min([A,B]))*10)/10;maxD=ceil(max(max([A,B]))*10)/10;
% else
%     minD=floor(min(min([B]))*10)/10;maxD=ceil(max(max([B]))*10)/10;
% end
% minD=min(minD,-maxD);maxD=max(maxD,-minD);
imagesc(B');
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize)
caxis([minC,maxC]);
title('$$B$$','interpreter','latex');
axis('square');
xlim([1,n-1]);
ylim([1,n-1]);
% clean_panel(ax,map2,pos,id,n,col,fontSize)


% C MCorr
% ax=subplot(s,t,2*t+3);
ax=subplot(s,t,2*t+2);
hold all
% MH=ceil(max(max(mcorrH(2:end,2:end))));
% mH=floor(min(min(mcorrH(2:end,2:end))));
% mH=min(mH,-MH);
% MH=max(MH,-mH);
imagesc(mcorrH');
caxis([mH,MH]);
set(gca,'YDir','normal')
title('$$C = A \circ B$$','FontSize',fontSize,'interpreter','latex');
axis('square');
xlim([1,n-1]);
ylim([1,n-1]);
% clean_panel(ax,map2,pos,id,n,col,fontSize)



%% MGC

warning('off','all');
% [~,~,~,optimalInd]=FindLargestRectangles((pMLocal<=pMGC), [0 0 1],[2,2]);
% optimalInd=find(optimalInd==1);
% if (pMLocal(end)<=pMGC && isempty(find(optimalInd==n*n, 1)))
%     optimalInd=n*n;
% end

% A MGC
% ax=subplot(s,t,4);
% ax=subplot('Position',[left(4), bottom(3), width, height]);
ax=subplot(s,t,3);
hold all
imagesc(A_MGC');
caxis([minC,maxC]);
set(gca,'YDir','normal')
title('$A^{k}$','interpreter','latex','FontSize',fontSize)
% title([{'3. MGC^k^,^l'}; {'(local scale)'}; {' '}],'FontSize',fontSize-1,...
%     'Units', 'normalized','Position', [0 1], 'HorizontalAlignment', 'left')
axis('square');
xlim([1,n-1]);
ylim([1,n-1]);
% clean_panel(ax,map2,pos,id,n,col,fontSize)


% B MGC
% ax=subplot(s,t,t+4);
ax=subplot(s,t,t+3);
hold all
imagesc(B_MGC');
caxis([minC,maxC]);
set(gca,'YDir','normal')
title('$$B^{l}$$','interpreter','latex');
axis('square');
xlim([1,n-1]);
ylim([1,n-1]);
% clean_panel(ax,map2,pos,id,n,col,fontSize)

% C MGC
% ax=subplot(s,t,2*t+4);
ax=subplot(s,t,2*t+3);
cla, hold all
imagesc(C_MGC');
set(gca,'YDir','normal')
caxis([mH,MH]);
title('$$C^{kl} = A^{k} \circ B^{l}$$','interpreter','latex');
axis('square');
xlim([1,n-1]);
ylim([1,n-1]);

ax=subplot(s,t,4);
hold all
plot(reshape(C,size(C,1)^2,1),reshape(D,size(D,1)^2,1),'k.','MarkerSize',3);
xlim([0,1]);
ylim([0,1]);
set(gca,'XTick',[0,1],'YTick',[0,1]);
title('Distance Scatter Plot','FontSize',fontSize, ...
    'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')
axis('square');

% B MGC
% ax=subplot(s,t,t+4);
map2 = brewermap(128,'BuPu'); % brewmap
ax=subplot(s,t,t+4);
hold all
map2=flipud(map2);
map2=flipud(map2);
colormap(ax,map2);
% cmap2=flipud(cmap2);
% set(groot,'defaultAxesColorOrder',cmap2);
% colormap(cmap2)
kmin=2;
ph=testMLocal';
%indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% ph(indPower)=2;
imagesc(ph);
plot(l,k,'go','markerSize',6,'linewidth',3)
plot(size(ph,2),size(ph,1),'.','markerSize',18,'MarkerFaceColor',glob,'Color',glob)
hold off
% colormap(ax,cmap2)
xlim([1 size(ph,2)]);
ylim([1 size(ph,1)]);
axis('square');
    title('Correlation Map','FontSize',fontSize, ...
        'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')

% C MGC
% ax=subplot(s,t,2*t+4);
ax=subplot(s,t,2*t+4);
hold on;
imagesc(log(pMLocal'));
    set(gca,'YDir','normal')
map2=flipud(map2);
colormap(ax,map2);
    set(gca,'FontSize',fontSize);
    %     if i==3
%     xlabel(xlabs{i},'FontSize',fs, ...
%        'Units', 'normalized','Position', [-0.01, -0.14], 'HorizontalAlignment', 'left')
%     ylabel(ylabs{i},'FontSize',fs, ...
%         'Units', 'normalized','Position', [-0.18 0], 'HorizontalAlignment', 'left')
    title('P-Value Map','FontSize',fontSize, ...
        'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')
    
    
    %[~,indP]=MGCScaleVerify(p2All',rep);
    indP=optimalInd;
    [m,n]=size(pMLocal);
    [J,I]=ind2sub([m,n],indP);
    Ymin=min(I);
    Ymax=max(I);
    Xmin=min(J);
    Xmax=max(J);
    %
    if Xmin==Xmax && Ymin==Ymax
         plot(Xmin,Ymin,'go','markerSize',6,'linewidth',3);
    else
        lw=2;
        plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
        plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
        plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
        plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
    end
    tmp=zeros(m,n);
    tmp(J,I)=1;
    tmp(testMLocal<testMGC)=0;
    [k,l]=ind2sub([m,n],find(tmp==1,1,'last'));
    plot(k,l,'go','markerSize',6,'linewidth',3);
    plot(m,n,'.','markerSize',18,'MarkerFaceColor',glob,'Color',glob)
    xticks=[5,round(m/2)-1,m-1];
%     if i==1,  xticks(1)=3; end
    %  set(gca,'XTick',xticks,'XTickLabel',[2,round(m/2),m]); % Remove x axis ticks
%    set(gca,'YTick',[3,round(n/2)-1,n-1],'YTickLabel',[2,round(n/2),n]); % Remove x axis ticks
    xlim([2,m]);
    ylim([2,n]);
    axis('square');
    hold off
        
    suptitle(fnames2{opt});
F.fname=pre2; %strcat(pre2, num2str(i));
F.wh=[11 9];
print_fig(gcf,F)

function [disRank]=DistRanks(dis)
% An auxiliary function that sorts the distance entries within each column by ascending order.
%
% The input is assumed to be a distance matrix.
%
% The output is column-wise rank, ordered from 1,...,n.
%
% For ties, the minimum ranking is used,
% e.g. if there are repeating distance entries, the order is like 1,2,3,3,4,..,n-1.

[n,m]=size(dis);
disRank=zeros(n,m);
for i=1:m
    [~,~,a]=unique(dis(:,i));
    disRank(:,i)=a;
end
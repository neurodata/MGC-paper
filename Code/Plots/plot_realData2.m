function []=plot_realData2(type)
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
pre2=strcat(rootDir,'Figures/FigReal');% The folder to save figures
sameBar=1;
if type==1
    load(strcat(rootDir,'Data/Preprocessed/BrainCP.mat'))
    load(strcat(rootDir,'Data/Results/CorrPermDistTestTypeBrainCxP.mat'))
    n=42;
    C=distC;
    D=distP;
end
%%%
%if select==2
%    load(strcat(rootDir,'Data/Preprocessed/BrainHippoShape.mat'))
%end
%%%
if type==2
    load(strcat(rootDir,'Data/Preprocessed/semipar.mat'))
    load(strcat(rootDir,'Data/Results/CorrPermDistTestTypeMigrainxCCI.mat'))
    n=109;
    distCCI=squareform(pdist(cci));
    C=distMigrain(ind,ind);
    D=distCCI(ind,ind);
end

if type==3
    load(strcat(rootDir,'Data/Preprocessed/proteomics.mat'))
%     per=(LabelIndAll~=2) & (LabelIndAll<5);
%     LabelIndAll(per)=1;
    per=(LabelIndAll<5);
    m=318;
    D=LabelIndAll(per);
    C=A(:,per)';
    D=squareform(pdist(D));
    C=squareform(pdist(C));
end

%% figure stuff

cmap=zeros(4,3);
gr =[0,1,0];
ma = [1,0,1];
map3 = brewermap(128,'PiYG'); % brewmap
lgr=map3(100,:);
dgr=map3(128,:);
gr=lgr;
glob= [0.5,0.5,0.5];
% cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = gr;
cmap(3,:) = gr;
cmap(4,:) = gr;
% cmap(3,:) = cy;
map1=cmap;
map2 = brewermap(128,'BuPu'); % brewmap
set(groot,'defaultAxesColorOrder',map2);

%% Mantel

if type~=3
    height=0.2;
    width=0.2;
    left=zeros(4,1);
    bottom=zeros(3,1);
    left(1)=0;
    left(2)=0.23;
    left(3)=0.46;
    left(4)=0.72;
    bottom(1)=0.7;
    bottom(2)=0.4;
    bottom(3)=0.1;
    fontSize=30;
%     s=3;t=4;
% A Mantel
% ax=subplot(s,t,2); %
[C,D,RX,RY]=MGCDistTransform(C,D,'mantel');
[m,n]=size(testMLocal);
ax=subplot('Position',[left(1), bottom(1), width, height]);
hold all
if sameBar==1
    maxC=ceil(max(max([C,D]))*10)/10;
    minC=ceil(min(min([C,D]))*10)/10;
else
    maxC=ceil(max(max(C))*10)/10;
    minC=ceil(min(min(C))*10)/10;
end
minC=min(minC,-maxC);maxC=max(maxC,-minC);
imagesc(C');
caxis([minC,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize); % Remove y axis ticks
%xlabel('# X Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
 %   'Position', [0 -0.16], 'HorizontalAlignment', 'left')
%ylabel('# Y Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
 %   'Position', [-0.25 0.4], 'VerticalAlignment', 'bottom')
% 'Position', [-0.2 -0.05], 'HorizontalAlignment', 'Left')
title('$$\tilde{A}$$','interpreter','latex');
% text(45,110,'$\tilde{A}$','interpreter','latex','FontSize',fontSize)
% title([{'1. Mantel'}; {'(pairwise distances)'}; {' ' }],'FontSize',fontSize-1, 'Units', 'normalized', ...
%     'Position', [0 1], 'HorizontalAlignment', 'left');
% title([{'Mantel'}; {' '}],'FontSize',fontSize);
%clean_panel(ax,map2,pos,id,n,col,fontSize)
set(gca,'visible','on')
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[1,round(m/2),m],'XtickLabel',[1,round(n/2),m]); % Remove x axis ticks
set(gca,'YTick',[1,round(n/2),n],'YtickLabel',[1,round(n/2),n]); % Remove x axis ticks
axis('square')
pos=get(ax,'position');
colorbar('location','eastoutside');
set(ax,'position',pos);


% B Mantel
% ax=subplot(s,t,t+2);
% ax=subplot(s,t,5);
ax=subplot('Position',[left(1), bottom(2), width, height]);
hold all
if sameBar~=1
    maxC=ceil(max(max(D))*10)/10;
end
imagesc(D');
set(gca,'FontSize',fontSize)
% colormap(ax,map3);
caxis([minC,maxC]);
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[],'YTick',[]); % Remove x axis ticks
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize); % Remove y axis ticks
title('$\tilde{B}$','interpreter','latex','FontSize',fontSize);
%clean_panel(ax,map2,pos,id,n,col,fontSize)
axis('square')

% C Mantel
% ax=subplot(s,t,2*t+2);
mantelH=C.*D;
% ax=subplot(s,t,9);
ax=subplot('Position',[left(1), bottom(3), width, height]);
hold all
MH=ceil(max(max(mantelH(2:end,2:end))));
mH=floor(min(min(mantelH(2:end,2:end))));
mH=min(mH,-MH);
MH=max(MH,-mH);
imagesc(mantelH');
axis('square')
% pos = get(ax,'position');
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[],'YTick',[]); % Remove x axis ticks
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize); % Remove y axis ticks
title('$$\tilde{C} = \tilde{A} \circ \tilde{B}$$','FontSize',fontSize,'interpreter','latex');
caxis([mH,MH]);
% colorbar('location','westoutside')
%clean_panel(ax,map2,pos,id,n,col,fontSize)
% set(ax,'position',pos);


%% Mcorr
[A,B]=MGCDistTransform(C,D,'mgc');
mcorrH=A.*B;
if sameBar==1
    minC=floor(min(min([A,B]))*10)/10;maxC=ceil(max(max([A,B]))*10)/10;
else
    minC=floor(min(min([A]))*10)/10;maxC=ceil(max(max([A]))*10)/10;
end
minC=min(minC,-maxC);maxC=max(maxC,-minC);


% A Mcorr
%  ax=subplot(s,t,2);
 ax=subplot('Position',[left(2), bottom(1), width, height]);
%ax=subplot('Position',[left(3), bottom(3), width, height]);
hold all
imagesc(A');
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize); % Remove y axis ticks
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[],'YTick',[]); % Remove x axis ticks
caxis([minC,maxC]);
title('$$A$$','interpreter','latex');
% text(45,110,'$A$','interpreter','latex','FontSize',fontSize)
% title([{'2. Mcorr'}; {'(single center)'}; {' '}],'FontSize',fontSize-1, 'Units', 'normalized', ...
%     'Position', [0 1], 'HorizontalAlignment', 'left')
%clean_panel(ax,map2,pos,id,n,col,fontSize)
axis('square')

% B MCorr
%  ax=subplot(s,t,6);
ax=subplot('Position',[left(2), bottom(2), width, height]);
%ax=subplot('Position',[left(3), bottom(2), width, height]);
hold all
if sameBar==1
    minD=floor(min(min([A,B]))*10)/10;maxD=ceil(max(max([A,B]))*10)/10;
else
    minD=floor(min(min([B]))*10)/10;maxD=ceil(max(max([B]))*10)/10;
end
minD=min(minD,-maxD);maxD=max(maxD,-minD);
imagesc(B');
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize)
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[],'YTick',[]); % Remove x axis ticks
caxis([minD,maxD]);
title('$$B$$','interpreter','latex');
%clean_panel(ax,map2,pos,id,n,col,fontSize)
axis('square')

% C MCorr
%  ax=subplot(s,t,10);
ax=subplot('Position',[left(2), bottom(3), width, height]);
%ax=subplot('Position',[left(3), bottom(1), width, height]);
hold all
MH=ceil(max(max(mcorrH(2:end,2:end)))*100)/100;
mH=floor(min(min(mcorrH(2:end,2:end)))*100)/100;
mH=min(mH,-MH);
MH=max(MH,-mH);
imagesc(mcorrH');
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[],'YTick',[]); % Remove x axis ticks
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize)
title('$$C = A \circ B$$','FontSize',fontSize,'interpreter','latex');
%clean_panel(ax,map2,pos,id,n,col,fontSize)
axis('square')


%% MGC

[k,l]=ind2sub(size(testMLocal),optimalInd(1));
% A MGC
%  ax=subplot(s,t,3);
ax=subplot('Position',[left(3), bottom(1), width, height]);
 %[A,B]=MGCDistTransform(C,D,'mgc');
 A_MGC=A;
 B_MGC=B;
 A_MGC(RX>k)=0;
 B_MGC(RY>l)=0;
 C_MGC=A_MGC.*B_MGC;
%ax=subplot('Position',[left(4), bottom(3), width, height]);
hold all
imagesc(A_MGC');
caxis([minC,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize)
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[],'YTick',[]); % Remove x axis ticks
title('$$A^{k}$$','interpreter','latex');
% text(45,110,'$A^{k}$','interpreter','latex','FontSize',fontSize)
% title([{'3. MGC^k^,^l'}; {'(local scale)'}; {' '}],'FontSize',fontSize-1,...
%     'Units', 'normalized','Position', [0 1], 'HorizontalAlignment', 'left')
%clean_panel(ax,map2,pos,id,n,col,fontSize)
axis('square')

% B MGC
% ax=subplot(s,t,7);
ax=subplot('Position',[left(3), bottom(2), width, height]);
%ax=subplot('Position',[left(4), bottom(2), width, height]);
hold all
imagesc(B_MGC');
caxis([minD,maxD]);
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize)
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[],'YTick',[]); % Remove x axis ticks
title('$$B^{l}$$','interpreter','latex');
%clean_panel(ax,map2,pos,id,n,col,fontSize)
axis('square')

% C MGC
% ax=subplot(s,t,11);
ax=subplot('Position',[left(3), bottom(3), width, height]);
%ax=subplot('Position',[left(4), bottom(1), width, height]);
cla, hold all
MH=ceil(max(max(C_MGC(2:end,2:end)))*100)/100;
mH=floor(min(min(C_MGC(2:end,2:end)))*100)/100;
mH=min(mH,-MH);MH=max(MH,-mH);
imagesc(C_MGC');
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize)
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[],'YTick',[]); % Remove x axis ticks
caxis([mH,MH]);
title('$$C^{kl} = A^{k} \circ B^{l}$$','interpreter','latex');
%clean_panel(ax,map2,pos,id,n,col,fontSize)
axis('square')

%% Col5 Multiscale Test Statistics
% ax=subplot(s,t,4);
ax=subplot('Position',[left(4), bottom(1), width, height]);
%ax=subplot('Position',[left(5), bottom(3), width, height]);
hold on
set(groot,'defaultAxesColorOrder',map1);
kmin=1;
ph=testMLocal(kmin:m,kmin:n)';
%indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% ph(indPower)=2;
imagesc(ph);
plot(m,n,'.','markerSize',30,'MarkerFaceColor',glob,'Color',glob)
plot(k,l,'go','markerSize',10,'linewidth',5)
hold off
axis('square')
set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
xlim([1 m]);
ylim([1 n]);
set(gca,'XTick',[1,round(m/2),m],'XtickLabel',[1,round(n/2),m]); % Remove x axis ticks
set(gca,'YTick',[1,round(n/2),n],'YtickLabel',[1,round(n/2),n]); % Remove x axis ticks
pos=get(ax,'position');
% cmap=map4;
colormap(ax,map2)
hm=ceil(max(max(ph))*100)/100;
hm=ceil(prctile(ph(ph<1),99)*100)/100;
caxis([0 hm])
h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
set(h,'FontSize',fontSize);
xlim([1 n]);
ylim([1 m]);
% set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(m/2),n],'YTickLabel',[2,round(m/2),n],'FontSize',fontSize);
% set(gca,'XTick',[],'YTick',[])
set(gca,'FontSize',fontSize)
set(ax,'position',pos);
% xlabel('# X Neighbors','FontSize',fontSize, ...
%     'Units', 'normalized','Position', [0 -0.2], 'HorizontalAlignment', 'left')
% ylabel('# Y Neighbors','FontSize',fontSize, ...
%     'Units', 'normalized','Position', [-0.2 0], 'HorizontalAlignment', 'left')
% text(-1,73,'4. Multiscale Maps','fontSize',fontSize,'fontweight','bold');
% title(1,60,[{'4. Multiscale Maps'}; {'(all scales)'}; {' '}],'FontSize',fontSize,...
%     'Units', 'normalized','Position', [0 1.1], 'HorizontalAlignment', 'left')
title('MGC Image');
% title('Local Correlations','fontweight','normal','FontSize',fontSize);
xlabel('# X Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
    'Position', [0 -0.2], 'HorizontalAlignment', 'left')
ylabel('# Y Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
    'Position', [-0.3 0.5], 'VerticalAlignment', 'bottom')
text(10,110,'Test Statistics','FontSize',fontSize)
%clean_panel(ax,map4,pos,id,n,col,fontSize)

% %% Col5 Multiscale Power Maps
% ax=subplot(s,t,12);
% ax=subplot('Position',[left(5), bottom(2), width, height]);
% hold on
% set(groot,'defaultAxesColorOrder',map1);
% kmin=2;
% ph=powerMLocal(kmin:n,kmin:n)';
% indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% % ph(indPower)=2;
% imagesc(ph);
% plot(n-1,n-1,'.','markerSize',30,'MarkerFaceColor',glob,'Color',glob)
% hold off
% 
% set(gca,'FontSize',fontSize)
% set(gca,'YDir','normal')
% cmap=map4;
% colormap(ax,cmap)
% caxis([0 1])
% h=colorbar('Ticks',[0,0.5,1]);%,'location','westoutside');
% set(h,'FontSize',fontSize);
% xlim([1 n-1]);
% ylim([1 n-1]);
% set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',16);
% set(gca,'XTick',[],'YTick',[])
% pos = get(ax,'position');
% %xlabel('# of X Neighbors','FontSize',fontSize, ...
% %      'Units', 'normalized','Position', [0 -0.2], 'HorizontalAlignment', 'left')
% %ylabel('# of Y Neighbors','FontSize',fontSize, ...
% %       'Units', 'normalized','Position', [-0.2 0], 'HorizontalAlignment', 'left')
% %text(-1,73,'4. Multiscale Maps','fontSize',fontSize,'fontweight','bold');
% %text(19,55,'Power','FontSize',fontSize)
% title('Powers','fontweight','normal','FontSize',fontSize);
% axis('square')
% clean_panel(ax,map4,pos,id,n,col,fontSize)

%% Col5 multiscale p-value map
% ax=subplot(s,t,12);
ax=subplot('Position',[left(4), bottom(3), width, height]);
%ax=subplot('Position',[left(5), bottom(1), width, height]);
hold on

%set(groot,'defaultAxesColorOrder',map1);
kmin=1;
% pMLocal(pMLocal<=eps)=0.005;
ph=pMLocal';
% ph(indP)=0.00001;
imagesc(log(ph)); %log(ph)-min(log(ph(:))));
pos = get(ax,'position');
set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize)
%cmap=map4;
colormap(ax,flipud(map2));
%ceil(max(max(ph))*10)/10
% caxis([0 1]);
cticks=[0.001, 0.01, 0.1, 0.5];

h=colorbar('Ticks',log(cticks),'TickLabels',cticks);%,'location','westoutside');
set(h,'FontSize',fontSize);
% set(gca,'XTick',[1,round(n/2)-1,n-1],'YTick',[1,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',16);
set(gca,'XTick',[],'YTick',[])
%xlabel('# of Neighbors for X','FontSize',16)
%ylabel('# of Neighbors for Y','FontSize',16) %,'Rotation',0,'position',[-7,20]);
% xlim([1 n-1]);
% ylim([1 n-1]);
% plot([n-2:n-2],[n-2:n-1],'-m','linewidth',2)
% plot([n-1:n-1],[n-2:n-1],'-m','linewidth',12)
% plot(k,l,'s','color',glob,'markerSize',5,'MarkerFaceColor',glob)
% plot(k,l,'s','color',loca,'markerSize',5,'MarkerFaceColor',loca)
plot(m,n,'.','markerSize',30,'MarkerFaceColor',glob,'Color',glob)
plot(k,l,'go','markerSize',10,'linewidth',5)

% draw boundary around optimal scale
%[pval,indP]=MGCScaleVerify(ph,1000);
indP=optimalInd(1);
%disp(strcat('Approximated MGC p-value: ',num2str(pval)));
% indP=indP(2:end,2:end)';
[J,I]=ind2sub(size(ph),indP);
Ymin=min(I);
Ymax=max(I);
Xmin=min(J);
Xmax=max(J);

lw=1.5;
plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
xlim([2,n]);
ylim([2,n]);
%     imagesc(k,l,1);
hold off
title('P-values');
% title('P-values','FontSize',fontSize, 'Units', 'normalized', ...
%     'Position', [0.2 1.2], 'HorizontalAlignment', 'center');
axis('square')
set(ax,'position',pos);
axis('square')

pre2=strcat(rootDir,'Figures/');% The folder to save figures
F.fname=strcat(pre2, strcat('RealDataFig',num2str(type)));
F.wh=[10 6.5]*2;
F.PaperPositionMode='auto';

print_fig(gcf,F)
end

if type==3;
    fontSize=15;
pos1 =[0,0.35,0.45,0.45];
ax=subplot('position',pos1);
hold on
kmin=1;
[testMGC, testMLocal, optimalInd]=MGCSampleStat(C,D);
[m,n]=size(testMLocal);
[k,l]=ind2sub(size(testMGC),optimalInd(1));
ph=testMLocal(kmin:m,kmin:n)';
%indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% ph(indPower)=2;
imagesc(ph);
plot(m,n,'.','markerSize',30,'MarkerFaceColor',glob,'Color',glob)
plot(l,k,'go','markerSize',10,'linewidth',5)
hold off
axis('square')
pos = get(ax,'position');

set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
% cmap=map4;
colormap(ax,map2)
hm=ceil(max(max(ph))*100)/100;
hm=ceil(prctile(ph(ph<1),99)*100)/100;
caxis([0 hm])
h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
set(h,'FontSize',fontSize);
% h2=get(h,'Position');
% h2(4)=h2(4)/2;
% set(h,'Position',h2)
xlim([1 m]);
ylim([1 n]);
%set(gca,'XTick',[2.5,round(n/2)-1,n],'YTick',[2.5,round(n/2)-1,n],'XTickLabel',[2,round(m/2),n],'YTickLabel',[2,round(m/2),n],'FontSize',16);
%set(gca,'XTick',[],'YTick',[])
% pos = get(ax,'position');
title('MGC Image','FontSize',fontSize-1, 'Units', 'normalized', ...
    'Position', [0 1.05], 'HorizontalAlignment', 'left');
% title('Local Correlations','fontweight','normal','FontSize',fontSize);
xlabel('# X Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
    'Position', [0 -0.17], 'HorizontalAlignment', 'left')
ylabel('# Y Neighbors','FontSize',fontSize, 'Units', 'normalized', ...
    'Position', [-0.15 0.3], 'VerticalAlignment', 'bottom')
%text(10,110,'Test Statistics','FontSize',fontSize)
% pos2 = get(ax,'position');
% pos2(3:4) = F.pos(3:4);
set(ax,'position',pos);

load(strcat(pre1,'ScreeningPancvsNormal.mat'))
testMGC2=testMGC;
load(strcat(pre1,'ScreeningPancvsAll.mat'))
pos1 =[0.52,0.35,0.45,0.45];
ax=subplot('position',pos1);
set(ax,'FontSize',fontSize-1);
% [~,ind]=sort(testMGC,'descend');
% plot(testMGC(ind),'.-','LineWidth',lw);
hold on
p1=testMGC2(:,6);
p2=testMGC(:,6);
plot(p1,p2,'.','Color',gr,'MarkerSize',15);
x=0:0.0001:1;
y=0.001;
plot(x,y*ones(length(x),1),'--','Color',glob);
plot(y*ones(length(x),1),x,'--','Color',glob);
ind=181;
plot(p1(ind),p2(ind),'.','Color',gr,'MarkerSize',25);
text(p1(ind),p2(ind),'neurogranin','VerticalAlignment','bottom','HorizontalAlignment','left','Color',gr,'FontSize',fontSize+5);
hold off
set(gca,'XScale','log','YScale','log');
% xlim([1,318]);
% ylim([-0.1,0.5]);
%set(gca,'YTick',[0,0.2,0.4],'XTick',[1,100,200,300]);
% ylabel('Magnitude','FontSize',fs, ...
%     'Units', 'normalized','Position', [-0.21 0], 'HorizontalAlignment', 'left')
% xlabel('Features','FontSize',fs,...
%     'Units', 'normalized','Position', [-0.01, -0.16], 'HorizontalAlignment', 'left')
ylabel('P-values for Panc vs All','FontSize',fontSize, ...
    'Units', 'normalized','Position', [-0.2 0], 'HorizontalAlignment', 'left')
xlabel('P-values for Panc vs Norm','FontSize',fontSize,...
    'Units', 'normalized','Position', [-0.01, -0.17], 'HorizontalAlignment', 'left')
title('Cancer Biomarker Discovery','FontSize',fontSize, ...
    'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left');
axis('square');
pre2=strcat(rootDir,'Figures/');% The folder to save figures
F.fname=strcat(pre2, strcat('RealDataFig',num2str(type)));
F.wh=[4 2.5]*2;
F.PaperPositionMode='auto';

print_fig(gcf,F)
end
%clean_panel(ax,map4,pos,id,n,col,fontSize)
%%
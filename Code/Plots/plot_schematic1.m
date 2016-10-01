function []=plot_schematic1(type)

% type=1;n=50;dim=1;noise=1;
% CorrSimPlotsA(type,n,dim,noise,pre1);
% Used to generate figure A in the files

%% % File path searching
if nargin<1
    type=1;
end

n=50;
noise=0;
if type==1
    noise=0.5;
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

figure('units','normalized','position',[0 0 1 1])
s=1;t=5;
try
    load(strcat(rootDir,'Data/Results/CorrFigure1Type',num2str(type),'n',num2str(n),'.mat')); % The folder to locate data
catch
    display('no file exist, running instead')
    run_fig1Data(type,n,noise);
    load(strcat(rootDir,'Data/Results/CorrFigure1Type',num2str(type),'n',num2str(n),'.mat')); % The folder to locate data
end

fontSize=18;
mkSize=20;
fontSize2=18;
tfs=18;

%% plotting parameters

cmap=zeros(2,3);
gray = [0.5,0.5,0.5];
black=[0,0,0];
map2 = brewermap(128,'PiYG'); % brewmap
map3 = map2(size(map2,1)/2+1:end,:);
% map3 = brewermap(128,'Greens'); % brewmap
map4 = brewermap(128,'GnBu'); % brewmap
gr=map2(120,:);
pu=map2(8,:);
loca=[0,1,0];
glob=[0.5,0.5,0.5];
mgc='Cyan';
mgc=loca;
cmap(1,:) = pu;
cmap(2,:) = gr;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

height=0.35; %18; %21;
width=0.16; %17;
hspace=0.04;
vspace=0.09;
for i=1:6
    left(i)=0.01+(i-1)*(width+hspace);
end
left(1)=0.025;
bottom=0.3;

% optimal scales
indP=optimalInd;
[J,I]=ind2sub(size(pMLocal),indP);
Ymin=min(I)-1;
Ymax=max(I)-1;
Xmin=min(J)-1;
Xmax=max(J)-1;
%%  Col 1
% ax=subplot(s,t,1);
ax=subplot('Position',[left(1), bottom, width, height]);
hold all
set(groot,'defaultAxesColorOrder',map2);
plot(x,y,'.','MarkerSize',mkSize,'Color',gray);
if type==1 && noise==1
xlabel('Cloud Shape','FontSize',fontSize,...
    'Units', 'normalized','Position', [-0.01, -0.2], 'HorizontalAlignment', 'left')
ylabel('Ground Wetness','FontSize',fontSize, ...
    'Units', 'normalized', 'Position', [-0.06 0], 'HorizontalAlignment', 'left')
else
xlabel('$x$','FontSize',fontSize2+6,'Interpreter','latex',...
    'Units', 'normalized','Position', [-0.01, -0.2], 'HorizontalAlignment', 'left')
ylabel('$y$','FontSize',fontSize2+6,'Interpreter','latex', ...
    'Units', 'normalized', 'Position', [-0.06 0], 'HorizontalAlignment', 'left')
end

if type>5
[I,J]=ind2sub([n,n],find(C_MGC>0.1,1,'first'));
J2=find(mcorrH(J,:)<0,1,'last');
else
    I=2;J=4;J2=n;
end
id=[I,J,J2,J];
id2=[1,2,3,2];
col=[0 0 0];
if abs(y(id(1))-y(id(2)))/(max(y)-min(y))<0.02
hy=[+5,-5,0]/100*(max(y)-min(y));
else
    hy=zeros(3,1);
end
hs=2/100*(max(x)-min(x));

for ind=[1,2,3]; %length(id)
    text(x(id(ind))+hs,y(id(ind))+hy(ind),num2str(ind),'fontsize',fontSize,'color',col)
    plot(x(id(ind)),y(id(ind)),'.','MarkerSize',mkSize,'Color',col);
end

tname=CorrSimuTitle(type);
findex=strfind(tname,'.');
tname=tname(findex+1:end);
xlim([min(x)-0.2, max(x)]);
ylim([min(y)-0.2, max(y)]);

if type == 1, AB='A'; else AB='B'; end
AB='';
tit1=strcat('0', AB ,'. Sample Data');
title([{tit1}; {' '}], 'Units', 'normalized', ...
    'Position', [0 1.1], 'HorizontalAlignment', 'left','FontSize',tfs);
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove x axis tick
axis('square')
pos = get(ax,'position');

%%  Pairwise Distances
ax=subplot('Position',[left(2), bottom, width, height]);
% ax=subplot(s,t,2);
x=reshape(C,n^2,1);
y=reshape(D,n^2,1);
RC=DistRanks(C);
RD=DistRanks(D)';
RC=(RC<=Xmax+1);
RD=(RD<=Ymax+1);
ind1=reshape(RC&RD,n^2,1);
hold on
set(groot,'defaultAxesColorOrder',map2);
plot(x(ind1==0),y(ind1==0),'.','MarkerSize',6,'Color',gray);
plot(x(ind1==1),y(ind1==1),'+','MarkerSize',4,'Color',loca);

x12=sub2ind([n,n], id(1),id(2));
x23=sub2ind([n,n], id(2),id(3));
text(x(x12)+hs,y(x12)+hy(1),'(1, 2)','fontsize',fontSize,'color',col)
plot(x(x12),y(x12),'.','MarkerSize',mkSize/2,'Color',col);

text(x(x23)+hs,y(x23)+hy(1),'(2, 3)','fontsize',fontSize,'color',col)
plot(x(x23),y(x23),'.','MarkerSize',mkSize/2,'Color',col);
hold off

tname=CorrSimuTitle(type);
findex=strfind(tname,'.');
tname=tname(findex+1:end);
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);
warning('off','all')
<<<<<<< HEAD
xlabel('${Distance}_x(x_i,x_j)$','FontSize',fontSize2+8,...
    'Units', 'normalized','Position', [-0.01, -0.06], 'HorizontalAlignment', 'left','Interpreter','latex');
ylabel('$Distance_y(y_i,y_j)$','FontSize',fontSize2+8, ...
=======
xlabel('Distance$$_{x}(x_i,x_j)$$','FontSize',fontSize2+4,...
    'Units', 'normalized','Position', [-0.01, -0.2], 'HorizontalAlignment', 'left','Interpreter','latex');
ylabel('Distance$$_{y}(y_i,y_j)$$','FontSize',fontSize2+4, ...
>>>>>>> 8b4c1eadc4187cf8228c1a7239f9eb6d21a2e6f1
    'Units', 'normalized', 'Position', [-0.06 0], 'HorizontalAlignment', 'left','Interpreter','latex');

tit1=strcat('1', AB ,'. Pairwise Distances');
title([{tit1}; {' '}], 'Units', 'normalized', ...
    'Position', [0 1.1], 'HorizontalAlignment', 'left','FontSize',tfs);
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove x axis tick
%pos=[nan, nan, width, height];
axis('square')
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

%% Multiscale Correlation Map %%%%%%%%%%%%
ax=subplot('Position',[left(3), bottom, width, height]);
% ax=subplot(s,t,3);
hold on
set(groot,'defaultAxesColorOrder',map1);
kmin=2;
ph=tA(kmin:n,kmin:n)';
%indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% ph(indPower)=2;
imagesc(ph);
[k,l]=find(tA==test);
plot(k-1,l-1,'gs','markerSize',5,'MarkerFaceColor','g')  
hold off

set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
cmap=map4;
colormap(cmap)
% hm=ceil(max(max(ph))*100)/100;
hm=ceil(prctile(ph(ph<1),99)*100)/100;
caxis([0 hm])
h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
set(h,'FontSize',fontSize);
axpos = ax.Position;
hpos=h.Position;
hpos(3)=0.5*hpos(3);
hpos(1)=hpos(1)+0.015;
h.Position=hpos;
ax.Position = axpos;
xlim([1 n-1]);
ylim([1 n-1]);
set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',fontSize);
xlabel('X Scales','FontSize',fontSize2,...
    'Units', 'normalized','Position', [-0.010, -0.20], 'HorizontalAlignment', 'left');
ylabel('Y Scales','FontSize',fontSize2, ...
    'Units', 'normalized', 'Position', [-0.22 -0.02], 'HorizontalAlignment', 'left');

tit1=strcat('2', AB ,'. Multiscale Correlation');
title([{tit1}; {'Map & Test Statistic'}],'FontSize',tfs, ...
   'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left','color','g')
axis('square')
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

%% Null Distributions
ax=subplot('Position',[left(4), bottom, width, height]);
% ax=subplot(s,t,4);
minp=min([min(tN(:,n,n)),min(tN(:,k,l)),tA(k,l),tN(n,n)]);
minp=floor(minp*10)/10;
maxp=max([max(tN(:,n,n)),max(tN(:,k,l)),tA(k,l),tN(n,n)]);
maxp=ceil(maxp*10)/10;
p=tN(:,k,l);
[f1,xi1]=ksdensity(p,'support',[-1,1]);
p=tN(:,n,n);
[f,xi]=ksdensity(p,'support',[-1,1]);
p=testN;
[f2,xi2]=ksdensity(p,'support',[-1,1]);

hold on
plot(xi,f,'.-','LineWidth',4,'Color',glob);
% plot(xi1,f1,'.-','LineWidth',4,'Color',loca);
plot(xi2,f2,'.-','LineWidth',4,'Color',mgc);
% legend('Mcorr','OMGC','SMGC','Location','East');
% legend boxoff
set(gca,'FontSize',fontSize);
x1=sum(sum(mcorrH))/norm(A,'fro')/norm(B,'fro');x2=sum(sum(C_MGC))/norm((A_MGC-mean(mean(A_MGC))),'fro')/norm((B_MGC--mean(mean(B_MGC))),'fro');
x1=round(x1*100)/100;x2=round(x2*100)/100;x3=round(test*100)/100;
plot(x1,0.1,'*','MarkerSize',12,'Color',glob,'linewidth',2);
% plot(x2,0.1,'*','MarkerSize',12,'Color',loca,'linewidth',2);
plot(x3,0.1,'*','MarkerSize',12,'Color',mgc,'linewidth',2);

% x1 = tA(end);
ind=find(xi>x1,1,'first');
y1=max(f)+2;
y2 = max(f1)+2;
y3 = 5;
txt1 = strcat('$$p(c) =', num2str(pMLocal(end)),'$$');
% txt2 = strcat('$$p(c^{*}) = ', num2str(pMLocal(k,l)),'$$');
txt3 = strcat('$$p(\hat{c}^{*}) = ', num2str(pMGC),'$$');
c=text(x3,y3,txt3,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',mgc,'Interpreter','latex');
set(c,'FontSize',fontSize);
% b=text(x2,y2,txt2,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',loca,'Interpreter','latex');
% set(b,'FontSize',fontSize2);
set(gca,'XTick',x3+0.1,'TickLength',[0 0],'XTickLabel',x3);
% if abs(x1-x3)>0.02 
a=text(x1,y1,txt1,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',glob,'Interpreter','latex');
set(a,'FontSize',fontSize);
 if abs(x1-x3)>0.02 
set(gca,'XTick',sort([x3+0.05,x1+0.05]),'TickLength',[0 0],'XTickLabel',sort([x3,x1]));
end
% if abs(x2-x3)>0.02 && abs(x2-x1)>0.02
% set(gca,'XTick',sort([x1+0.02,x2+0.02,x3+0.02]),'TickLength',[0 0],'XTickLabel',sort([x1,x2,x3]));
% end
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',fontSize);
ylim([0 y1+10]);

xlim([minp,min(maxp+0.1,1)]);
if type<5
    xlim([minp,min(maxp+0.4,1)]);
end
xlabel('Test Statistic','FontSize',fontSize2,...
    'Units', 'normalized','Position', [-0.01, -0.20], 'HorizontalAlignment', 'left');
ylabel('Density','FontSize',fontSize2, ...
    'Units', 'normalized', 'Position', [-0.10, 0], 'HorizontalAlignment', 'left');
set(gca,'YTick',[])

tit1=strcat('3', AB ,'. Null Distributions');
title([{tit1}; {'& P-Values'}],'FontSize',tfs, ...
   'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left')
axis('square')
hold off
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

%% Multiscale P-Value Map
ax=subplot('Position',[left(5), bottom, width, height]);
% ax=subplot(s,t,5);
hold on
set(groot,'defaultAxesColorOrder',map1);
kmin=2;
ph=pMLocal(kmin:n,kmin:n)';
imagesc(log(ph));

% draw boundary around optimal scale
% indP=optimalInd;
% [J,I]=ind2sub(size(pMLocal),indP);
% Ymin=min(I)-1;
% Ymax=max(I)-1;
% Xmin=min(J)-1;
% Xmax=max(J)-1;
lw=3;
plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
xlim([2,n]);
ylim([2,n]);
%     imagesc(k,l,1);
hold off

set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
cmap=map4;
colormap(ax,flipud(cmap));
%ceil(max(max(ph))*10)/10
% caxis([0 1]);
cticks=[0.001, 0.01, 0.1, 0.5];
h=colorbar('Ticks',log(cticks),'TickLabels',cticks);%,'location','eastoutside');
set(h,'FontSize',fontSize);
axpos = ax.Position;
hpos=h.Position;
hpos(3)=0.5*hpos(3);
hpos(1)=hpos(1)+0.023;
h.Position=hpos;
ax.Position = axpos;
xlim([1 n-1]);
ylim([1 n-1]);
set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',fontSize);
xlabel('X Scales','FontSize',fontSize2,...
    'Units', 'normalized','Position', [-0.010, -0.20], 'HorizontalAlignment', 'left');
ylabel('Y Scales','FontSize',fontSize2, ...
    'Units', 'normalized', 'Position', [-0.22 -0.02], 'HorizontalAlignment', 'left');

tit1=strcat('4', AB ,'. Multiscale P-Value');
title([{tit1}; {'Map & Optimal Scales'}],'FontSize',tfs, ...
   'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left','color','g')
axis('square')
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
h=suptitle(CorrSimuTitle(type));
set(h,'FontSize',fontSize2+10,'Units', 'normalized','Position', [0.01, -0.60,0], 'HorizontalAlignment', 'left','rotation',90)

%
pre2=strcat(rootDir,'Figures/');% The folder to save figures
donzo=0;
if donzo==1
    F.fname=strcat(pre2, 'Fig',num2str(type));
else
    F.fname=strcat(pre2, 'Auxiliary/A2_type', num2str(type),'_n', num2str(n), '_noise', num2str(round(noise*10)));
end
F.wh=[10 3]*2;
F.PaperPositionMode='auto';

print_fig(gcf,F)
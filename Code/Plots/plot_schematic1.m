function []=plot_schematic1(type)

% type=1;n=50;dim=1;noise=1;
% CorrSimPlotsA(type,n,dim,noise,pre1);
% Used to generate figure A in the files

%% % File path searching
if nargin<1
    type=1;
    n=50;
    noise=1;
else
    type=8;
    n=50;
    noise=0;
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

figure('units','normalized','position',[0 0 1 1])
s=1;t=3;
try
    load(strcat(rootDir,'Data/Results/CorrFigure1Type',num2str(type),'n',num2str(n),'.mat')); % The folder to locate data
catch
    display('no file exist, running instead')
    run_fig1Data(type,n,noise);
    load(strcat(rootDir,'Data/Results/CorrFigure1Type',num2str(type),'n',num2str(n),'.mat')); % The folder to locate data
end

fontSize=18;
mkSize=20;
sameBar=0;


%% plotting parameters

cmap=zeros(2,3);
gray = [0.5,0.5,0.5];
map2 = brewermap(128,'PiYG'); % brewmap
map3 = map2(size(map2,1)/2+1:end,:);
% map3 = brewermap(128,'Greens'); % brewmap
map4 = brewermap(128,'GnBu'); % brewmap
gr=map2(120,:);
pu=map2(8,:);
loca=[0,1,0];
glob=[0.5,0.5,0.5];
mgc='Cyan';
cmap(1,:) = pu;
cmap(2,:) = gr;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

% % [left,bottom,width,height]
% height=0.21; %18; %21;
% width=height-0.04; %17;
% hspace=0;
% vspace=0.09;
% for i=1:6
%     left(i)=0.01+(i-1)*(width+hspace);
%     bottom(i)=0.06+(i-1)*(height+vspace);
% end
% bottom(3)=bottom(3)+0.01;
% left(5)=left(5)+0.03;
% left(6)=left(6)+0.09;
% left(2:end)=left(2:end)+0.02;
% 
% fig=figure(1); clf
% set(gcf,'units','normalized','position',[0 0 1 1])


%%  Col 1
% ax=subplot(s,t,t+1);
ax=subplot(s,t,1);
hold all
set(groot,'defaultAxesColorOrder',map2);
plot(x,y,'.','MarkerSize',mkSize,'Color',gray);
if type==1 && noise==1
xlabel('Cloud Shape','FontSize',fontSize,...
    'Units', 'normalized','Position', [-0.008, -0.1], 'HorizontalAlignment', 'left')
ylabel('Ground Wetness','FontSize',fontSize, ...
    'Units', 'normalized', 'Position', [-0.03 0], 'HorizontalAlignment', 'left')
else
xlabel('X','FontSize',fontSize,...
    'Units', 'normalized','Position', [-0.008, -0.1], 'HorizontalAlignment', 'left')
ylabel('Y','FontSize',fontSize, ...
    'Units', 'normalized', 'Position', [-0.03 0], 'HorizontalAlignment', 'left')
end

if type>5
[I,J]=ind2sub([n,n],find(C_MGC>0.1,1,'first'));
J2=find(mcorrH(J,:)<0,1,'last');
else
    I=1;J=5;J2=20;
end
id=[I,J,J2,J];
id2=[1,2,3,2];
col=[1 .5 0];
hy=[-5,+5,0]/100*(max(y)-min(y));
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

% title(['0. ', [tname], ' (X,Y)'], 'Units', 'normalized', ...
title([{'Sample Data'}], 'Units', 'normalized', ...
    'Position', [0 1.1], 'HorizontalAlignment', 'left')
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove x axis tick
%pos=[nan, nan, width, height];
axis('square')
pos = get(ax,'position');

ax=subplot(s,t,2);
%ax=subplot('Position',[left(1), bottom(1)+width/2+0.01, width, height]);

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
plot(xi1,f1,'.-','LineWidth',4,'Color',loca);
plot(xi2,f2,'.-','LineWidth',4,'Color',mgc);
set(gca,'FontSize',fontSize);
% x1=round(tA(end)*100)/100;
x1=sum(sum(mcorrH))/norm(A,'fro')/norm(B,'fro');x2=sum(sum(C_MGC))/norm((A_MGC-mean(mean(A_MGC))),'fro')/norm((B_MGC--mean(mean(B_MGC))),'fro');
x1=round(x1*100)/100;x2=round(x2*100)/100;x3=round(test*100)/100;
plot(x1,0.1,'*','MarkerSize',12,'Color',glob,'linewidth',2);
plot(x2,0.1,'*','MarkerSize',12,'Color',loca,'linewidth',2);
plot(x3,0.1,'*','MarkerSize',12,'Color',mgc,'linewidth',2);
h=legend('Mcorr','OMGC','SMGC','Location','NorthEast');
legend boxoff
set(h,'FontSize',fontSize-3);
% set(gca,'XTickLabel',[x1;x2],'YTick',[]); % Remove x axis ticks

% x1 = tA(end);
ind=find(xi>x1,1,'first');
% x2 = tA(k,l);
% x3=test;
y1=max(f)+2;
y2 = max(f1)+2;
y3 = 5;
%txt1 = {'Mcorr';['p = ' num2str(pMLocal(end))]};
txt1 = strcat('$$p(c) =', num2str(pMLocal(end)),'$$');
txt2 = strcat('$$p(c^{*}) = ', num2str(pMLocal(k,l)),'$$');
txt3 = strcat('$$p(\hat{c}^{*}) = ', num2str(pMGC),'$$');
c=text(x3,y3,txt3,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',mgc,'Interpreter','latex');
set(c,'FontSize',fontSize);
set(gca,'XTick',x3,'TickLength',[0 0],'XTickLabel',x3);
if abs(x2-x3)>0.03
b=text(x2,y2,txt2,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',loca,'Interpreter','latex');
set(b,'FontSize',fontSize);
set(gca,'XTick',sort([x3+0.02,x2+0.02]),'TickLength',[0 0],'XTickLabel',sort([x3,x2]));
end
if abs(x1-x3)>0.03 && abs(x1-x2)>0.03
a=text(x1,y1,txt1,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',glob,'Interpreter','latex');
set(a,'FontSize',fontSize);
set(gca,'XTick',sort([x1+0.02,x2+0.02,x3+0.02]),'TickLength',[0 0],'XTickLabel',sort([x1,x2,x3]));
end
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',fontSize-6);
ylim([0 y1+10]);
% <<<<<<< HEAD
xlim([minp,maxp+0.1]);
xlabel('Test Statistic','FontSize',fontSize,...
    'Units', 'normalized','Position', [-0.008, -0.1], 'HorizontalAlignment', 'left')
ylabel('Density','FontSize',fontSize, ...
    'Units', 'normalized', 'Position', [-0.03 0], 'HorizontalAlignment', 'left')
% =======
% set(a,'FontSize',fontSize);
% set(b,'FontSize',fontSize);
% set(c,'FontSize',fontSize);
% xlim([minp,maxp+0.1]);
% xlabel('Test Statistic','FontSize',fontSize-5,'HorizontalAlignment','right');
% ylabel('Density','FontSize',fontSize-5, ...
% >>>>>>> c55135b820b0d669c0dcdac7dc915b0200706c61
%     'Units', 'normalized', 'Position', [-0.02 0], 'HorizontalAlignment', 'left')
set(gca,'YTick',[])
title([{'Tests by Mcorr and its MGC'}],'FontSize',fontSize, ...
   'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left')
% pos2 = get(ax,'position');
% pos2(3:4) = [pos(3:4)];
% set(ax,'position',pos2);
axis('square')
hold off
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

ax=subplot(s,t,3);
hold on
set(groot,'defaultAxesColorOrder',map1);
kmin=2;
ph=tA(kmin:n,kmin:n)';
%indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% ph(indPower)=2;
imagesc(ph);

% draw boundary around optimal scale
%[pval,indP]=MGCScaleVerify(ph,1000);
indP=optimalInd;
%disp(strcat('Approximated MGC p-value: ',num2str(pval)));
% indP=indP(2:end,2:end)';
[J,I]=ind2sub(size(ph),indP);
Ymin=min(I);
Ymax=max(I);
Xmin=min(J);
Xmax=max(J);

plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',3)
plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',3)
plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',3)
plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',3)
xlim([2,n]);
ylim([2,n]);
%     imagesc(k,l,1);
hold off

set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
cmap=map4;
colormap(cmap)
%hm=ceil(max(max(ph))*100)/100;
hm=ceil(prctile(ph(ph<1),99)*100)/100;
caxis([0 hm])
h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
set(h,'FontSize',fontSize-6);
xlim([1 n-1]);
ylim([1 n-1]);
set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',fontSize-6);
xlabel('# of X Neighbors','FontSize',fontSize,...
    'Units', 'normalized','Position', [-0.008, -0.1], 'HorizontalAlignment', 'left')
ylabel('# of Y Neighbors','FontSize',fontSize, ...
    'Units', 'normalized', 'Position', [-0.1 0], 'HorizontalAlignment', 'left')
title([{'Multiscale Correlation Map'}],'FontSize',fontSize, ...
   'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left')
% xlabel('# of X Neighbors','FontSize',fontSize, ...
%     'Units', 'normalized','Position', [0 -0.2], 'HorizontalAlignment', 'left')
% ylabel('# of Y Neighbors','FontSize',fontSize, ...
%     'Units', 'normalized','Position', [-0.2 0], 'HorizontalAlignment', 'left')
% text(-1,73,'4. Multiscale Maps','fontSize',fontSize,'fontweight','bold');
% title(1,60,[{'4. Multiscale Maps'}; {'(all scales)'}; {' '}],'FontSize',fontSize,...
%     'Units', 'normalized','Position', [0 1.1], 'HorizontalAlignment', 'left')
%text(-1,70,[{'4. Multiscale Maps'};{'(all scales)'}],'fontSize',fontSize,'fontweight','bold');
% title('Local Correlations','fontweight','normal','FontSize',fontSize);
%text(10,55,'Test Statistics','FontSize',fontSize)
axis('square')
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

%
pre2=strcat(rootDir,'Figures/');% The folder to save figures
donzo=1;
if donzo==1
    F.fname=strcat(pre2, 'Fig',num2str(type));
else
    F.fname=strcat(pre2, 'Auxiliary/A2_type', num2str(type),'_n', num2str(n), '_noise', num2str(round(noise*10)));
end
F.wh=[8 3]*2;
F.PaperPositionMode='auto';

print_fig(gcf,F)
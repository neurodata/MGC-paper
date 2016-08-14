function []=plot_schematic(newSim,type,n,noise)

% type=6;n=100;dim=1;noise=1;
% CorrSimPlotsA(type,n,dim,noise,pre1);
% Used to generate figure A in the files

%% % File path searching
if nargin<1
    newSim=0;
end
if nargin<2
    type=13;
end
if nargin<3
    n=50;
end
if nargin<4
    noise=0;
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

if newSim==1
    run_fig1Data(type,n,noise);
end
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
s=3;t=5;
gr=map2(120,:);
pu=map2(8,:);
loca=[0,1,0];
glob= [1,0,1];
cmap(1,:) = pu;
cmap(2,:) = gr;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

% [left,bottom,width,height]
height=0.21; %18; %21;
width=height-0.04; %17;
hspace=0;
vspace=0.09;
for i=1:6
    left(i)=0.01+(i-1)*(width+hspace);
    bottom(i)=0.06+(i-1)*(height+vspace);
end
bottom(3)=bottom(3)+0.01;
left(5)=left(5)+0.03;
left(6)=left(6)+0.09;
left(2:end)=left(2:end)+0.02;

figure(1), clf
set(gcf,'units','normalized','position',[0 0 1 1])

%% permute sample order

xynorm=sum([x,y]').^2';
[foo, bar]=sort(xynorm);
bar=1:n;
%bar=randperm(n);


%%  Col 1
% ax=subplot(s,t,t+1); 
ax=subplot('Position',[left(1), bottom(2)+width/2+0.01, width, height]);
cla, hold all
set(groot,'defaultAxesColorOrder',map2);
plot(x,y,'.','MarkerSize',mkSize,'Color',gray);
xlabel('x')
ylabel('y')

[I,J]=ind2sub([n,n],find(C_MGC>0.1,1,'first'));
J2=find(mcorrH(J,:)<0,1,'last');
ids=unique([I,J]);
% xx=15;
id=[I,J,J2,J];
id2=[1,2,3,2];
col=[1 .5 0];
hy=[-0.4,0.4,0];

for ind=[1,2,3]; %length(id)
    hs=0.2;
    text(x(id(ind))+hs,y(id(ind))+hy(ind),num2str(ind),'fontsize',fontSize,'color',col)
    plot(x(id(ind)),y(id(ind)),'.','MarkerSize',mkSize,'Color',col);
end

tname=CorrSimuTitle(type);
findex=strfind(tname,'.');
tname=tname(findex+1:end);

% title(['0. ', [tname], ' (X,Y)'], 'Units', 'normalized', ...
title([{'0. Sample Data'}], 'Units', 'normalized', ...
'Position', [0 1.1], 'HorizontalAlignment', 'left')
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove x axis tick
pos=[nan, nan, width, height];
axis('square')

% make table
clc

MantelVec=[C(id(1),id(2)), D(id(1),id(2)), mantelH(id(1),id(2)), C(id(3),id(4)), D(id(3),id(4)), mantelH(id(3),id(4)), mantelH(id(1),id(2))+mantelH(id(3),id(4))];
McorrVec=[A(id(1),id(2)), B(id(1),id(2)), mcorrH(id(1),id(2)), A(id(3),id(4)), B(id(3),id(4)), mcorrH(id(3),id(4)), mcorrH(id(1),id(2))+mcorrH(id(3),id(4))];
MGCVec=[A_MGC(id(1),id(2)), B_MGC(id(1),id(2)), A_MGC(id(1),id(2))*B_MGC(id(1),id(2)), A_MGC(id(3),id(4)), B_MGC(id(3),id(4)), A_MGC(id(3),id(4)) * B_MGC(id(3),id(4)), A_MGC(id(1),id(2))*B_MGC(id(1),id(2))+A_MGC(id(3),id(4)) * B_MGC(id(3),id(4))];

[MantelVec; McorrVec; MGCVec]'
display(['k=', num2str(k), ' l=',num2str(l)])

formatSpec = 'k= %1i, l= %1i\n\n';
fprintf(formatSpec,k,l) 

formatSpec = '\n & Mantel & Mcorr & MGC \\\\ \n';
fprintf(formatSpec)
formatSpec = '\n $\\delta_x$(%i,%i) & %1.2f & %1.2f & %1.2f  \\\\ \n $\\delta_y$(%i,%i) & %1.2f & %1.2f & %1.2f  \\\\ \n $\\delta_x \\times \\delta_y$ & %1.2f & %1.2f & %1.2f  \\\\ \n \n';
fprintf(formatSpec,id2(1),id2(2),C(id(1),id(2)),A(id(1),id(2)),A_MGC(id(1),id(2)),...
                   id2(1),id2(2),D(id(1),id(2)),B(id(1),id(2)),B_MGC(id(1),id(2)),...
                   mantelH(id(1),id(2)),mcorrH(id(1),id(2)),A_MGC(id(1),id(2))*B_MGC(id(1),id(2)))
fprintf('\\hline\n\n')
formatSpec = '\n $\\delta_x$(%i,%i) & %1.2f & %1.2f & %1.2f  \\\\ \n $\\delta_y$(%i,%i) & %1.2f & %1.2f & %1.2f  \\\\ \n $\\delta_x \\times \\delta_y$ & %1.2f & %1.2f & %1.2f  \\\\ \n \n';
fprintf(formatSpec,id2(3),id2(4),C(id(3),id(4)),A(id(3),id(4)),A_MGC(id(3),id(4)),...
                   id2(3),id2(4),D(id(3),id(4)),B(id(3),id(4)),B_MGC(id(3),id(4)),...
                   mantelH(id(3),id(4)),mcorrH(id(3),id(4)),A_MGC(id(3),id(4)) * B_MGC(id(3),id(4)))

fprintf('\\hline\n\n')

formatSpec = '\n $\\sum \\delta_x \\times \\delta_y$ & %1.2f & %1.2f & %1.2f  \\\\ \n \n';
fprintf(formatSpec,sum(sum(mantelH)), sum(sum(mcorrH)),sum(sum(C_MGC)));
fprintf(formatSpec,sum(sum(mantelH))/norm(C,'fro')/norm(D,'fro'), sum(sum(mcorrH))/norm(A,'fro')/norm(B,'fro'),sum(sum(C_MGC))/norm((A_MGC-mean(mean(A_MGC))),'fro')/norm((B_MGC--mean(mean(B_MGC))),'fro'));


%% Mantel


% A Mantel
% ax=subplot(s,t,2); %  
ax=subplot('Position',[left(2), bottom(3), width, height]);
hold all
if sameBar==1
    maxC=ceil(max(max([C,D]))*10)/10;
    minC=ceil(min(min([C,D]))*10)/10;
else
    maxC=ceil(max(max(C))*10)/10;
    minC=ceil(min(min(C))*10)/10;
end
imagesc(C(bar,bar)');
caxis([minC,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',16); % Remove y axis ticks
xlabel('sample index','FontSize',fontSize, 'Units', 'normalized', ...
'Position', [0 -0.15], 'HorizontalAlignment', 'left')
ylabel('sample index','FontSize',fontSize, 'Units', 'normalized', ...
'Position', [-0.3 0.72], 'VerticalAlignment', 'Top')
% 'Position', [-0.2 -0.05], 'HorizontalAlignment', 'Left')
text(24,53,'$\tilde{A}$','interpreter','latex','FontSize',fontSize)
title([{'1. Mantel'}; {'(pairwise dist''s)'}; {' ' }],'FontSize',fontSize, 'Units', 'normalized', ...
'Position', [0 1.1], 'HorizontalAlignment', 'left');
% title([{'Mantel'}; {' '}],'FontSize',fontSize);
clean_panel(ax,map2,pos,id,n,col,fontSize)
set(gca,'visible','on')
set(gca,'XTick',[1,round(n/2),n]); % Remove x axis ticks
set(gca,'YTick',[1,round(n/2),n]); % Remove x axis ticks

% B Mantel
% ax=subplot(s,t,t+2); 
ax=subplot('Position',[left(2), bottom(2), width, height]);
hold all
if sameBar~=1
    maxC=ceil(max(max(D))*10)/10;
end
imagesc(D(bar,bar)');
set(gca,'FontSize',16)
% colormap(ax,map3);
caxis([minC,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',13); % Remove y axis ticks
title('$\tilde{B}$','interpreter','latex','FontSize',fontSize);
clean_panel(ax,map2,pos,id,n,col,fontSize)


% C Mantel
% ax=subplot(s,t,2*t+2); 
ax=subplot('Position',[left(2), bottom(1), width, height]);
hold all
MH=ceil(max(max(mantelH(2:end,2:end))));
mH=floor(min(min(mantelH(2:end,2:end))));
mH=min(mH,-MH);
MH=max(MH,-mH);
imagesc(mantelH(bar,bar)');
set(gca,'YDir','normal')
title('$$\tilde{A} \circ \tilde{B}$$','FontSize',fontSize,'interpreter','latex');
caxis([mH,MH]);
colorbar('location','westoutside')
clean_panel(ax,map2,pos,id,n,col,fontSize)


%% Mcorr
if sameBar==1
    minC=floor(min(min([A,B]))*10)/10;maxC=ceil(max(max([A,B]))*10)/10;
else
    minC=floor(min(min([A]))*10)/10;maxC=ceil(max(max([A]))*10)/10;
end
minC=min(minC,-maxC);maxC=max(maxC,-minC);


% A Mcorr
% ax=subplot(s,t,3); 
ax=subplot('Position',[left(3), bottom(3), width, height]);
hold all
imagesc(A(bar,bar)');
set(gca,'YDir','normal')
caxis([minC,maxC]);
text(24,53,'$A$','interpreter','latex','FontSize',fontSize)
title([{'2. Mcorr'}; {'(single center)'}; {' '}],'FontSize',fontSize, 'Units', 'normalized', ...
'Position', [0 1.1], 'HorizontalAlignment', 'left')
clean_panel(ax,map2,pos,id,n,col,fontSize)


% B MCorr
% ax=subplot(s,t,t+3); 
ax=subplot('Position',[left(3), bottom(2), width, height]);
hold all
if sameBar==1
    minD=floor(min(min([A,B]))*10)/10;maxD=ceil(max(max([A,B]))*10)/10;
else
    minD=floor(min(min([B]))*10)/10;maxD=ceil(max(max([B]))*10)/10;
end
minD=min(minD,-maxD);maxD=max(maxD,-minD);
imagesc(B(bar,bar)');
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize)
caxis([minD,maxD]);
title('$$B$$','interpreter','latex');
clean_panel(ax,map2,pos,id,n,col,fontSize)


% C MCorr
% ax=subplot(s,t,2*t+3); 
ax=subplot('Position',[left(3), bottom(1), width, height]);
hold all
MH=ceil(max(max(mcorrH(2:end,2:end))));
mH=floor(min(min(mcorrH(2:end,2:end))));
mH=min(mH,-MH);
MH=max(MH,-mH);
imagesc(mcorrH(bar,bar)');
set(gca,'YDir','normal')
title('$$A \circ B$$','FontSize',fontSize,'interpreter','latex');
clean_panel(ax,map2,pos,id,n,col,fontSize)



%% MGC


% A MGC
% ax=subplot(s,t,4);
ax=subplot('Position',[left(4), bottom(3), width, height]);
hold all
imagesc(A_MGC(bar,bar)');
caxis([minC,maxC]);
set(gca,'YDir','normal')
text(24,53,'$A$','interpreter','latex','FontSize',fontSize)
title([{'3. MGC^k^,^l'}; {'(local scales)'}; {' '}],'FontSize',fontSize,...
    'Units', 'normalized','Position', [0 1.1], 'HorizontalAlignment', 'left')
clean_panel(ax,map2,pos,id,n,col,fontSize)


% B MGC
% ax=subplot(s,t,t+4); 
ax=subplot('Position',[left(4), bottom(2), width, height]);
hold all
imagesc(B_MGC(bar,bar)');
caxis([minD,maxD]);
set(gca,'YDir','normal')
title('$$B^{l^{*}}$$','interpreter','latex');
clean_panel(ax,map2,pos,id,n,col,fontSize)

% C MGC
% ax=subplot(s,t,2*t+4); 
ax=subplot('Position',[left(4), bottom(1), width, height]);
cla, hold all
MH=ceil(max(max(C_MGC(2:end,2:end))));
mH=floor(min(min(C_MGC(2:end,2:end))));
mH=min(mH,-MH);MH=max(MH,-mH);
imagesc(C_MGC(bar,bar)');
set(gca,'YDir','normal')
caxis([mH,MH]);
title('$$A^{k^{*}} \circ B^{l^{*}}$$','interpreter','latex');
clean_panel(ax,map2,pos,id,n,col,fontSize)



%% Col5 Multiscale Power Maps
% ax=subplot(s,t,2*t+t);
ax=subplot('Position',[left(5), bottom(3), width, height]);
hold on
set(groot,'defaultAxesColorOrder',map1);
kmin=2;
ph=power1(kmin:n,kmin:n)';
imagesc(ph);
hold off

set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
cmap=map4;
colormap(ax,cmap)
caxis([0 1])
h=colorbar('Ticks',[0,0.5,1]);%,'location','westoutside');
set(h,'FontSize',fontSize);
xlim([1 n-1]);
ylim([1 n-1]);
set(gca,'XTick',[1,round(n/2)-1,n-1],'YTick',[1,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',16);
pos = get(ax,'position');
xlabel('# of Neighbors for X','FontSize',fontSize)
ylabel('# of Neighbors for Y','FontSize',fontSize) %,'Rotation',0,'position',[-7,20]);
text(-1,73,'4. Multiscale Maps','fontSize',fontSize,'fontweight','bold');
text(19,55,'Power','FontSize',fontSize)
axis('square')


%% Col5 multiscale p-value map
% ax=subplot(s,t,t+t);
ax=subplot('Position',[left(5), bottom(2), width, height]);
hold on

set(groot,'defaultAxesColorOrder',map1);
kmin=2;
pAll(pAll<=eps)=0.005;
ph=pAll(2:end,2:end)';
% ph(ph<=eps)=0.0005;
imagesc(log(ph)); %log(ph)-min(log(ph(:))));

set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
cmap=map4;
colormap(ax,flipud(cmap));
%ceil(max(max(ph))*10)/10
% caxis([0 1]);
cticks=[0.001, 0.01, 0.1, 0.5];
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

h=colorbar('Ticks',log(cticks),'TickLabels',cticks);%,'location','westoutside');
set(h,'FontSize',fontSize);
% set(gca,'XTick',[1,round(n/2)-1,n-1],'YTick',[1,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',16);
set(gca,'XTick',[],'YTick',[])
%xlabel('# of Neighbors for X','FontSize',16)
%ylabel('# of Neighbors for Y','FontSize',16) %,'Rotation',0,'position',[-7,20]);
xlim([1 n-1]);
ylim([1 n-1]);
% plot([n-2:n-2],[n-2:n-1],'-m','linewidth',2)
% plot([n-1:n-1],[n-2:n-1],'-m','linewidth',12)
plot(n-1,n-1,'.m','markersize',24)


%     imagesc(k,l,1);
hold off
% pos2 = get(ax,'position');
% pos2(3:4) = [pos(3:4)];
% set(ax,'position',pos2);
%     %set(h,'FontSize',16,'location','southoutside');
%set(gca,'XTick',[]); % Remove x axis ticks
%set(gca,'YTick',[]); % Remove y axis ticks
title('P-Values','fontweight','normal','FontSize',fontSize);
%colorbar('location','westoutside');
axis('square')

%% Col 5 p-value
% subplot(s,t,t)
ax=subplot('Position',[left(5), bottom(1), width, height]);

minp=min([min(tN(:,n,n)),min(tN(:,k,l)),tA(k,l),tN(n,n)]);
minp=floor(minp*10)/10;
maxp=max([max(tN(:,n,n)),max(tN(:,k,l)),tA(k,l),tN(n,n)]);
maxp=ceil(maxp*10)/10;
p=tN(:,k,l);
[f1,xi1]=ksdensity(p,'support',[-1,1]);
p=tN(:,n,n);
[f,xi]=ksdensity(p,'support',[-1,1]);

hold on
plot(xi,f,'.-','LineWidth',4,'Color',glob);
plot(xi1,f1,'.-','LineWidth',4,'Color',loca);
set(gca,'FontSize',15);
% x1=round(tA(end)*100)/100;
x1=sum(sum(mcorrH))/norm(A,'fro')/norm(B,'fro');x2=sum(sum(C_MGC))/norm((A_MGC-mean(mean(A_MGC))),'fro')/norm((B_MGC--mean(mean(B_MGC))),'fro');
x1=round(x1*100)/100;x2=round(x2*100)/100;
plot(x1,0.1,'*','MarkerSize',12,'Color',glob,'linewidth',2);
plot(x2,0.1,'*','MarkerSize',12,'Color',loca,'linewidth',2);
set(gca,'XTick',[x1+0.05,x2+0.05],'TickLength',[0 0],'XTickLabel',[x1,x2]);
% set(gca,'XTickLabel',[x1;x2],'YTick',[]); % Remove x axis ticks

x1 = tA(end);
ind=find(xi>x1,1,'first');
x2 = tA(k,l);
y1=max(f)+2;
y2 = max(f1)+2;
txt1 = {'Mcorr';['p = ' num2str(pAll(end))]};
if pAll(k,l)<0.001;
    txt2 = {'MGC';'p < 0.001'};
else
    txt2 = {'MGC';['p = ' num2str(pAll(k,l))]};
end
a=text(x1,y1,txt1,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',glob);
b=text(x2,y2,txt2,'VerticalAlignment','middle','HorizontalAlignment','left','Color',loca);
ylim([0 y1+25]);
set(a,'FontSize',fontSize);
set(b,'FontSize',fontSize);
xlim([minp,maxp+0.04]);
xlabel('Test Statistic','FontSize',fontSize,'HorizontalAlignment','right');
ylabel('Density','FontSize',fontSize, ...
    'Units', 'normalized', 'Position', [-0.02 0], 'HorizontalAlignment', 'left')
set(gca,'YTick',[])
title([{'5. Test using'}; {'Optimal Scales'}],'FontSize',fontSize, ...
    'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left')
% pos2 = get(ax,'position');
% pos2(3:4) = [pos(3:4)];
% set(ax,'position',pos2);
% axis('square')
hold off
% axis([0, y2, 0, y2])
% 1-pAll(end)
% 1-pAll(k,l)
% power1(end)
% power1(k,l)


%%
pre2=strcat(rootDir,'Figures/');% The folder to save figures
donzo=1;
if donzo==1
    F.fname=strcat(pre2, 'FigA');    
else
    F.fname=strcat(pre2, 'Auxiliary/A2_type', num2str(type),'_n', num2str(n), '_noise', num2str(round(noise*10)));
end
F.wh=[10 6.5]*2;F.fname

print_fig(gcf,F)
function []=plot_auxiliary(newSim,type,n,noise)

% type=6;n=100;dim=1;noise=1;
% CorrSimPlotsA(type,n,dim,noise,pre1);
% Used to generate figure A in the files

%%% File path searching
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
end
pre2=strcat(rootDir,'Figures/Fig');% The folder to save figures

cc=1;
fontSize=20;
mkSize=20;
sameBar=0;

cmap=zeros(2,3);
gray = [0.5,0.5,0.5];
% ma = [1,0,1];
% % cy = [0,1,1];
% cmap(1,:) = ma;
% cmap(2,:) = ma;
% cmap(3,:) = ma;
% cmap(4,:) = ma;
% cmap(3,:) = cy;
% map0=gr;
% map1=cmap;
map2 = brewermap(128,'PRGn'); % brewmap
map3 = brewermap(128,'Greens'); % brewmap
map4 = brewermap(128,'GnBu'); % brewmap
s=3;t=4;
gr=map2(100,:);
pu=map2(28,:);
cmap(1,:) = pu;
cmap(2,:) = gr;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

figure(1), clf
set(gcf,'units','normalized','position',[0 0 1 1])

%%  Col 1
ax=subplot(s,t,1);
set(groot,'defaultAxesColorOrder',map2);
plot(x,y,'.','MarkerSize',mkSize,'Color',gray);
% xlim([-1.2,1.2]);
% ylim([-1.2,1.2]);
% axis('tight')

title('Scatter Plot of (X,Y)')
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove x axis tick
pos = get(ax,'position');
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')

%% Col 2
ax=subplot(s,t,2);
if sameBar==1
    maxC=ceil(max(max([C,D]))*10)/10;
else
    maxC=ceil(max(max(C))*10)/10;
end
imagesc(C');
colormap(ax,map3);
caxis([0,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',13); % Remove y axis ticks
xlabel('i','FontSize',fontSize,'position',[-10,0]);
ylabel('j','Rotation',0,'position',[-25,40],'FontSize',fontSize);
set(gca,'XTick',[round(n/2)-1,n-1],'XTickLabel',[round(n/2),n]); % Remove x axis ticks
set(gca,'YTick',[round(n/2)-1,n-1],'YTickLabel',[round(n/2),n]); % Remove x axis ticks
% title('Mantel $\tilde{A}=\delta_{x}(x_{i},x_{j})$','interpreter','latex','FontSize',fontSize);
title('Mantel','FontSize',fontSize);
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')

ax=subplot(s,t,3);
if sameBar~=1
    maxC=ceil(max(max(D))*10)/10;
end
imagesc(D');
set(gca,'FontSize',16)
colormap(ax,map3);
caxis([0,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',13); % Remove y axis ticks
%set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
xlabel('i','FontSize',fontSize,'position',[-10,0]);
ylabel('j','Rotation',0,'position',[-25,40],'FontSize',fontSize);
set(gca,'XTick',[round(n/2)-1,n-1],'XTickLabel',[round(n/2),n]); % Remove x axis ticks
set(gca,'YTick',[round(n/2)-1,n-1],'YTickLabel',[round(n/2),n]); % Remove x axis ticks
title('$\tilde{B}=\delta_{y}(y_{i},y_{j})$','interpreter','latex','FontSize',fontSize);
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')

%% Centering
RC=disToRanks(C);
RD=disToRanks(D);
if option==3;
    A=sum(sum(C))/n/(n-1);
    B=sum(sum(D))/n/(n-1);
    A=C-A;
    B=D-B;
    for j=1:n
        A(j,j)=0;
        B(j,j)=0;
    end
else
    H=eye(n)-(1/n)*ones(n,n);
    A=H*C;
    B=D*H;
    % % For mcorr, further adjust the double centered matrices to remove high-dimensional bias
    if option==2
        A=A-C/n;
        B=B-D/n;
        for j=1:n
            A(j,j)=0;
            B(j,j)=0;
        end
    end
end
%global stat
if cc==1
    A=A-mean(mean(A));B=B-mean(mean(B));
end
mcorrH=A.*B;
%

%% Col 3
ax=subplot(s,t,6);
if sameBar==1
    minC=floor(min(min([A,B]))*10)/10;maxC=ceil(max(max([A,B]))*10)/10;
else
    minC=floor(min(min([A]))*10)/10;maxC=ceil(max(max([A]))*10)/10;
end
minC=min(minC,-maxC);maxC=max(maxC,-minC);
imagesc(A');
set(gca,'YDir','normal')
colormap(ax,map2)
caxis([minC,maxC]);
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize);
% title('Mcorr $$A$$','FontSize',fontSize,'interpreter','latex');
title('Mcorr','FontSize',fontSize);
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')

ax=subplot(s,t,7);
if sameBar==1
    minD=floor(min(min([A,B]))*10)/10;maxD=ceil(max(max([A,B]))*10)/10;
else
    minD=floor(min(min([B]))*10)/10;maxD=ceil(max(max([B]))*10)/10;
end
minD=min(minD,-maxD);maxD=max(maxD,-minD);
imagesc(B');
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize)
colormap(ax,map2)
caxis([minD,maxD]);
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
title('$$B$$','interpreter','latex');
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')

%% Global to Local
RC=(RC>k);
RD=(RD>l);
A(RC)=0;
B(RD)=0;
if cc==1
    A=A-mean(mean(A));B=B-mean(mean(B));
end


%% Col 4
ax=subplot(s,t,5);
imagesc(A');
colormap(ax,map2)
caxis([minC,maxC]);
set(gca,'YDir','normal')
% title('MGC $$A^{k^{*}}$$','interpreter','latex');
title('MGC');
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')

ax=subplot(s,t,10);
imagesc(B');
set(gca,'FontSize',fontSize)
colormap(ax,map2)
caxis([minD,maxD]);
set(gca,'YDir','normal')
title('$$B^{l^{*}}$$','interpreter','latex');
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')

%% Col3 Center
ax=subplot(s,t,4);
MH=ceil(max(max(mcorrH(2:end,2:end)))*10)/10;
mH=floor(min(min(mcorrH(2:end,2:end)))*10)/10;
mH=min(mH,-MH);MH=max(MH,-mH);
imagesc(mcorrH');
set(gca,'YDir','normal')
colormap(ax,map2)
colorbar('Ticks',[mH,0,MH]);
caxis([mH,MH]);
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
title('$$A \circ B$$','FontSize',fontSize,'interpreter','latex');
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')

%% %Col4 Center
ax=subplot(s,t,11);
mcorrH=A.*B;
MH=ceil(max(max(mcorrH(2:end,2:end)))*10)/10;
mH=floor(min(min(mcorrH(2:end,2:end)))*10)/10;
mH=min(mH,-MH);MH=max(MH,-mH);
imagesc(mcorrH');
set(gca,'YDir','normal')
colormap(ax,map2)
caxis([mH,MH]);
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
title('$$A^{k^{*}} \circ B^{l^{*}}$$','interpreter','latex');
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')

%% Col5 MPM
set(groot,'defaultAxesColorOrder',map1);
ax=subplot(s,t,12);
kmin=2;
hold on
ph=power1(kmin:n,kmin:n)';
imagesc(ph);
set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
cmap=map4;
colormap(ax,cmap)
caxis([0 1])
h=colorbar('Ticks',[0,0.5,1]);%,'location','westoutside');
set(h,'FontSize',fontSize);
set(gca,'XTick',[round(n/2)-1,n-1],'YTick',[round(n/2)-1,n-1],'XTickLabel',[round(n/2),n],'YTickLabel',[round(n/2),n],'FontSize',16);
set(gca,'XTick',[]);
xlabel('# of Neighbors for X','FontSize',16)
ylabel('# of Neighbors for Y','FontSize',16) %,'Rotation',0,'position',[-7,20]);
xlim([1 n-1]);
ylim([1 n-1]);

pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

%     imagesc(k,l,1);
hold off
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
%     %set(h,'FontSize',16,'location','southoutside');
%set(gca,'XTick',[]); % Remove x axis ticks
%set(gca,'YTick',[]); % Remove y axis ticks
title('Multiscale Power Map','FontSize',16);
axis('square')


%% Col5 p-value map
set(groot,'defaultAxesColorOrder',map1);
ax=subplot(s,t,8);
kmin=2;
hold on
ph=pAll(2:end,2:end)';
%ph(ph<eps)=eps;
imagesc(ph); %log(ph)-min(log(ph(:))));
set(gca,'FontSize',fontSize)
set(gca,'YDir','normal')
cmap=map4;
colormap(ax,flipud(cmap));
%ceil(max(max(ph))*10)/10
caxis([0 0.2]);
h=colorbar('Ticks',[0 0.05 0.2]);%,'location','westoutside');
set(h,'FontSize',fontSize);
set(gca,'XTick',[round(n/2)-1,n-1],'YTick',[round(n/2)-1,n-1],'XTickLabel',[round(n/2),n],'YTickLabel',[round(n/2),n],'FontSize',16);
set(gca,'XTick',[]);
%xlabel('# of Neighbors for X','FontSize',16)
%ylabel('# of Neighbors for Y','FontSize',16) %,'Rotation',0,'position',[-7,20]);
xlim([1 n-1]);
ylim([1 n-1]);

pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

%     imagesc(k,l,1);
hold off
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
%     %set(h,'FontSize',16,'location','southoutside');
%set(gca,'XTick',[]); % Remove x axis ticks
%set(gca,'YTick',[]); % Remove y axis ticks
title('Multiscale P-Value Map','FontSize',16);
%colorbar('location','westoutside');
axis('square')

%% Col 5 p-value
subplot(s,t,9)
minp=min([min(tN(:,n,n)),min(tN(:,k,l)),tA(k,l),tN(n,n)]);
minp=floor(minp*10)/10;
maxp=max([max(tN(:,n,n)),max(tN(:,k,l)),tA(k,l),tN(n,n)]);
maxp=ceil(maxp*10)/10;
p=tN(:,k,l);
[f1,xi1]=ksdensity(p,'support',[-1,1]);
p=tN(:,n,n);
[f,xi]=ksdensity(p,'support',[-1,1]);

hold on
plot(xi,f,'.-',xi1,f1,'.-','LineWidth',4);
set(gca,'FontSize',15);
x1=round(tA(end)*100)/100;
x2=round(tA(k,l)*100)/100;
plot(x1,0.1,'.','MarkerSize',mkSize,'Color',pu);
plot(x2,0.1,'*','MarkerSize',10,'Color',gr);
set(gca,'XTick',[x1+0.04,x2+0.04],'TickLength',[0 0]);
set(gca,'XTickLabel',[x1;x2],'YTick',[]); % Remove x axis ticks

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
a=text(x1,y1,txt1,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',pu);
b=text(x2,y2,txt2,'VerticalAlignment','middle','HorizontalAlignment','left','Color',gr);
xm=min([xi1,xi]);
xM=max([xi1,xi]);
xm=min([x1,x2,xm]);
xM=max([x1,x2,xM]);
xlim([xm xM]);
ylim([0 y1+25]);
set(a,'FontSize',17);
set(b,'FontSize',17);
xlim([minp,maxp]);
xlabel('Test Statistic','FontSize',fontSize);
ylabel('Density','FontSize',fontSize);
title('Test Statistics Distributions','FontSize',16);
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
axis('square')
hold off
% axis([0, y2, 0, y2])
% 1-pAll(end)
% 1-pAll(k,l)
% power1(end)
% power1(k,l)
k
l

%%
F.fname=strcat(pre2, 'A2_type', num2str(type),'_n', num2str(n), '_noise', num2str(round(noise*10)));
F.wh=[8 5]*2;
print_fig(gcf,F)
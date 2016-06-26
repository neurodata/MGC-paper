function []=plot_auxiliary

% type=6;n=100;dim=1;noise=1;
% CorrSimPlotsA(type,n,dim,noise,pre1);
% Used to generate figure A in the files

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));
load(strcat(rootDir,'Data/Results/CorrFigure1.mat')); % The folder to locate data
pre2=strcat(rootDir,'Figures/Fig');% The folder to save figures

cc=1;
fontSize=20;
mkSize=20;

cmap=zeros(4,3);
gr = [0.5,0.5,0.5];
ma = [1,0,1];
% cy = [0,1,1];
cmap(1,:) = ma;
cmap(2,:) = ma;
cmap(3,:) = ma;
cmap(4,:) = ma;
% cmap(3,:) = cy;
map0=gr;
map1=cmap;
map2 = brewermap(128,'PRGn'); % brewmap
map3 = brewermap(128,'Greens'); % brewmap
map4 = brewermap(128,'GnBu'); % brewmap
set(groot,'defaultAxesColorOrder',map1);
s=3;t=4;

figure('units','normalized','position',[0 0 1 1])
%%% Col 1
ax=subplot(s,t,1);
set(groot,'defaultAxesColorOrder',map0);
    plot(x,y,'.','MarkerSize',mkSize);
xlim([-1.2,1.2]);
ylim([-1.2,1.2]);

title('Scatter Plot of (X,Y)')
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove x axis tick
pos = get(ax,'position');
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

%%% Col 2
ax=subplot(s,t,2);
    maxC=ceil(max(max([C,D]))*10)/10;
imagesc(C');
colormap(ax,map3);
caxis([0,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',13); % Remove y axis ticks
xlabel('i','FontSize',fontSize,'position',[-5,0]);
ylabel('j','Rotation',0,'position',[-7,25],'FontSize',fontSize);
title('$\tilde{A}=\delta_{x}(x_{i},x_{j})$','interpreter','latex','FontSize',fontSize);
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

ax=subplot(s,t,3);
imagesc(D');
set(gca,'FontSize',16)
colormap(ax,map3);
caxis([0,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',13); % Remove y axis ticks
%set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
xlabel('i','FontSize',fontSize,'position',[-5,0]);
ylabel('j','Rotation',0,'position',[-7,25],'FontSize',fontSize);
title('$\tilde{B}=\delta_{y}(y_{i},y_{j})$','interpreter','latex','FontSize',fontSize);
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

%Centering
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

%%% Col 3
ax=subplot(s,t,6);
minC=floor(min(min([A,B]))*10)/10;maxC=ceil(max(max([A,B]))*10)/10;
minC=min(minC,-maxC);maxC=max(maxC,-minC);
imagesc(A');
set(gca,'YDir','normal')
colormap(ax,map2)
caxis([minC,maxC]);
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize);
title('$$A$$','FontSize',fontSize,'interpreter','latex');
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

ax=subplot(s,t,7);
minD=floor(min(min([A,B]))*10)/10;maxD=ceil(max(max([A,B]))*10)/10;
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

% Global to Local
RC=(RC>k);
RD=(RD>l);
A(RC)=0;
B(RD)=0;
if cc==1
    A=A-mean(mean(A));B=B-mean(mean(B));
end

%%% Col 4
ax=subplot(s,t,9);
imagesc(A');
colormap(ax,map2)
caxis([minC,maxC]);
set(gca,'YDir','normal')
title('$$A^{k^{*}}$$','interpreter','latex');
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

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

%Col3 Center
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
title('$$A.*B$$','FontSize',fontSize,'interpreter','latex');
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);

%Col4 Center
ax=subplot(s,t,8);
mcorrH=A.*B;
imagesc(mcorrH');
set(gca,'YDir','normal')
colormap(ax,map2)
caxis([mH,MH]);
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
title('$$A^{k^{*}}.*B^{l^{*}}$$','interpreter','latex');
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
% Col5
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
h=colorbar('Ticks',[0,0.5,1]);
set(h,'FontSize',fontSize);
set(gca,'XTick',[10,20,30,40],'YTick',[10,20,30,40],'FontSize',16);
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
colorbar

% Col 5 p-value
subplot(s,t,5)
minp=min([min(tN(:,n,n)),min(tN(:,k,l)),tA(k,l),tN(n,n)]);
minp=floor(minp*10)/10;
maxp=max([max(tN(:,n,n)),max(tN(:,k,l)),tA(k,l),tN(n,n)]);
maxp=ceil(maxp*10)/10;
p=tN(:,k,l);
[f1,xi1]=ksdensity(p,'support',[-1,1]);
p=tN(:,n,n);
[f,xi]=ksdensity(p,'support',[-1,1]);

hold on
plot(xi,f,'.:',xi1,f1,'.-','LineWidth',4);
set(gca,'FontSize',15);
plot(tA(end),0.1,'.','MarkerSize',mkSize);
plot(tA(k,l),0.1,'*','MarkerSize',10);
set(gca,'XTick',[round(tA(end)*100)/100,round(tA(k,l)*100)/100]);
set(gca,'XTickLabel',['      0.03'; '      0.31'],'YTick',[]); % Remove x axis ticks

x1 = tA(end);
ind=find(xi>x1,1,'first');
x2 = tA(k,l);
y1=max(f);
y2 = max(f1);
txt1 = {'Mcorr';['p = ' num2str(1-pAll(end))]};
txt2 = {'MGC';'p < 0.001'};
a=text(x1,y1,txt1,'VerticalAlignment','middle','HorizontalAlignment','left');
b=text(x2,y2,txt2,'VerticalAlignment','middle','HorizontalAlignment','left');
ylim([0 35]);
set(a,'FontSize',17);
set(b,'FontSize',17);
xlim([minp,maxp]);
ylabel('Density','FontSize',fontSize);
title('Test Statistics Distributions','FontSize',16);
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
hold off
1-pAll(end)
1-pAll(k,l)
power1(end)
power1(k,l)

%
F.fname=strcat(pre2, 'A2');
F.wh=[8 5]*2;
print_fig(gcf,F)
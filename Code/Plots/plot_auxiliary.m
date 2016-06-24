function []=plot_auxiliary(type)

% type=6;n=100;dim=1;noise=1;
% CorrSimPlotsA(type,n,dim,noise,pre1);
% Used to generate figure A in the files

%%%
fpath = mfilename('fullpath');
findex=strfind(fpath,'\');
rootDir=fpath(1:findex(end-2));
p = genpath(rootDir);
gits=strfind(p,'.git');
colons=strfind(p,':');
for i=0:length(gits)-1
    endGit=find(colons>gits(end-i),1);
    p(colons(endGit-1):colons(endGit)-1)=[];
end
addpath(p);

if nargin<1
    type=11;
end
pre1='..\..\Data\Results\'; % The folder to locate data
pre2='..\..\Figures\Fig'; % The folder to save figures
option=2;
n=50;
dim=1;
noise=0;
cc=1;
rep=1000;
repp=1;
fontSize=20;
mkSize=30;

cmap=zeros(3,3);
gr = [0.5,0.5,0.5];
ma = [1,0,1];
% cy = [0,1,1];
cmap(1,:) = ma;
cmap(2,:) = ma;
% cmap(3,:) = cy;
map1=cmap;
map2 = brewermap(128,'PRGn'); % brewmap
map3 = brewermap(128,'Greens'); % brewmap
map4 = brewermap(128,'GnBu'); % brewmap
set(groot,'defaultAxesColorOrder',map1);
s=3;t=4;

% Generate data
[x, y]=CorrSampleGenerator(type,n,dim,1, noise);
if noise~=0
    [x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
end
% Get distance matrix
[x,ind]=sort(x,'ascend');
y=y(ind);
C=squareform(pdist(x));
D=squareform(pdist(y));

% optimal scale
tA=zeros(n,n,rep);
tN=zeros(n,n,rep);
pa=zeros(n,n);
for rr=1:repp
    for r=1:rep;
        [x, y]=CorrSampleGenerator(type,n,dim,1, noise);
        CA=squareform(pdist(x));
        DA=squareform(pdist(y));
        tA(:,:,r)=LocalCorr(CA,DA,2);
        [x, y]=CorrSampleGenerator(type,n,dim,0, noise);
        CA=squareform(pdist(x));
        DA=squareform(pdist(y));
        tN(:,:,r)=LocalCorr(CA,DA,2);
    end
    power1=zeros(n,n);
    alpha=0.05;
    for i=1:n;
        for j=1:n;
            dCorT=sort(tN(i,j,:),'descend');
            cut1=dCorT(ceil(rep*alpha));
            power1(i,j)=mean(tA(i,j,:)>cut1);
        end
    end
    pa=pa+power1/repp;
end
power1=pa;
power1(1,:)=0;power1(:,1)=0; % Set the powers of all local tests at rank 0 to 0
neighbor=maxNeighbors(power1,0,tN,tA);

% Permutation p-value
tA=LocalCorr(C,D,2);
tN=zeros(rep,n,n);
pAll=zeros(n,n);
for r=1:rep;
    per=randperm(n);
    tmp=LocalCorr(C,D(per,per),2);
    tN(r,:,:)=tmp;
    if r==1
        pAll=(tmp<tA)/rep;
    else
        pAll=pAll+(tmp<tA)/rep;
    end
end

l=ceil(neighbor/n)
k=neighbor-n*(l-1)
minp=min([min(tN(:,n,n)),min(tN(:,k,l)),tA(k,l),tN(n,n)]);
minp=floor(minp*10)/10;
maxp=max([max(tN(:,n,n)),max(tN(:,k,l)),tA(k,l),tN(n,n)]);
maxp=ceil(maxp*10)/10;
p=tN(:,k,l);
[f1,xi1]=ksdensity(p,'support',[-1,1]);
p=tN(:,n,n);
[f,xi]=ksdensity(p,'support',[-1,1]);


figure('units','normalized','position',[0 0 1 1])
%%% Col 1
subplot(s,t,1);
if noise==0
    plot(x,y,'.','MarkerSize',mkSize);
else
    plot(x,y,'.',x1,y1,'k.','MarkerSize',mkSize);
end
%xlabel('X');
%ylabel('Y','Rotation',0,'position',[-1.4,-0.1]);

    xlim([-1.2,1.2]);
    ylim([-1.2,1.2]);

title('Scatter Plot of (X,Y)')
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove x axis tick

%%% Col 2
ax=subplot(s,t,2);
    maxC=ceil(max(max([C,D]))*10)/10;
imagesc(C');
colormap(ax,map3);
caxis([0,maxC]);
set(gca,'YDir','normal')
set(gca,'FontSize',fontSize); % Remove y axis ticks
xlabel('i')
ylabel('j')
title('$\tilde{A}=\delta_{x}(x_{i},x_{j})$','interpreter','latex');
pos = get(ax,'position');

ax=subplot(s,t,3);
imagesc(D');
set(gca,'FontSize',16)
colormap(ax,map3);
caxis([0,maxC]);
set(gca,'YDir','normal')
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
xlabel('i')
ylabel('j')
title('$\tilde{B}=\delta_{y}(y_{i},y_{j})$','interpreter','latex');

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

ax=subplot(s,t,10);
imagesc(B');
set(gca,'FontSize',fontSize)
colormap(ax,map2)
caxis([minD,maxD]);
set(gca,'YDir','normal')
title('$$B^{l^{*}}$$','interpreter','latex');
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks

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
% pos2 = get(ax,'position');
% pos2(3:4) = [pos(3:4)];
% set(ax,'position',pos2);

%Col4 Center
ax=subplot(s,t,8);
mcorrH=A.*B;
imagesc(mcorrH');
set(gca,'YDir','normal')
colormap(ax,map2)
caxis([mH,MH]);
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
title('$$A^{k^{*}}.*B^{l^{*}}$$','interpreter','latex');

% Col5
ax=subplot(s,t,12);
kmin=2;
hold on
ph=power1(kmin:n,kmin:n)';
imagesc(ph);
set(gca,'FontSize',16)
set(gca,'YDir','normal')
cmap=map4;
colormap(ax,cmap)
caxis([0.1 1])
xlabel('Number Of Neighbors for X','FontSize',13)
ylabel('Number Of Neighbors for Y','FontSize',13) %,'Rotation',0,'position',[-7,20]);
xlim([1 n-1]);
ylim([1 n-1]);
h=colorbar('Ticks',[0.1,1]);
set(h,'FontSize',fontSize);

% pos2 = get(ax,'position');
% pos2(3:4) = [pos(3:4)];
% set(ax,'position',pos2);

%     imagesc(k,l,1);
hold off
%     %set(h,'FontSize',16,'location','southoutside');
%set(gca,'XTick',[]); % Remove x axis ticks
%set(gca,'YTick',[]); % Remove y axis ticks
title('Multiscale Power Map','FontSize',fontSize);
colorbar

% Col 5 p-value
subplot(s,t,5)
hold on
plot(xi,f,'.:',xi1,f1,'.-','LineWidth',3);
set(gca,'FontSize',15);
plot(tA(end),0.1,'.','MarkerSize',mkSize);
plot(tA(k,l),0.1,'*','MarkerSize',10);
set(gca,'XTick',[0,tA(k,l)],'YTick',[]); % Remove x axis ticks

x1 = tA(end);
y1 = 0.01;
x2 = tA(k,l);
y2 = 0.01;
txt1 = strcat('Mcorr p = ',num2str(1-pAll(end)));
txt2 = strcat('MGC p < 0.001');
a=text(x1,y1,txt1,'VerticalAlignment','top','HorizontalAlignment','center');
b=text(x2,y2,txt2,'VerticalAlignment','top','HorizontalAlignment','center');
set(a,'FontSize',13);
set(b,'FontSize',13);
xlim([minp,maxp]);
ylabel('Probability Density','FontSize',15);
title('Test Statistics Distributions','FontSize',fontSize);

hold off
1-pAll(end)
1-pAll(k,l)
power1(end)
power1(k,l)

%
F.fname=strcat(pre2, 'A2');
F.wh=[8 4]*2;
print_fig(gcf,F)
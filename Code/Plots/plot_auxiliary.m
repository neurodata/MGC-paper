function []=plot_auxiliary(type)%,pre2)

% type=6;n=100;dim=1;noise=1;pre1='../../Data/';
% CorrSimPlotsA(type,n,dim,noise,pre1);
% Used to generate figure A in the files

if nargin<1
    type=11;
end
pre1='../../Data/Results/'; % The folder to locate data
pre2='../../Figures/Fig'; % The folder to save figures
% if nargin<2
%     pre2='../../Draft/Figures/FigReal'; % The folder to save figures
% end
option=2;
n=50;
dim=1;
noise=0;
cc=1;
samebar=1;
colorb=0;
rep=1000;
repp=1;
fontSize=20;
mkSize=20;
cmap=zeros(3,3);
gr = [0.5,0.5,0.5];
ma = [1,0,1];
% cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = gr;
% cmap(3,:) = cy;
map1=cmap;
map2 = brewermap(128,'PRGn'); % brewmap
map3 = brewermap(128,'Greens'); % brewmap
map4 = brewermap(128,'GnBu'); % brewmap
set(groot,'defaultAxesColorOrder',map1);
s=3;t=4;

figure('units','normalized','position',[0 0 1 1])
subplot(s,t,1);
% Generate data
[x, y]=CorrSampleGenerator(type,n,dim,1, noise);
if noise~=0
    [x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
end
% A1: scatter plot
if noise==0
    plot(x,y,'.','MarkerSize',mkSize);
else
    plot(x,y,'.',x1,y1,'k.','MarkerSize',mkSize);
end
xlabel('X');
ylabel('Y','Rotation',0,'position',[-1.4,-0.1]);
if samebar==1
    xlim([-1.2,1.2]);
    ylim([-1.2,1.2]);
end
title('Scatter Plot of (X,Y)')
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove x axis ticks
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
% figure
% ksdensity(reshape(tA(n,n,:),1,rep));
% figure
% ksdensity(reshape(tN(n,n,:),1,rep));
power1(1,:)=0;power1(:,1)=0; % Set the powers of all local tests at rank 0 to 0

neighbor=maxNeighbors(power1,0,tN,tA);
% [~,pAll]=MGCPermutationTest(C,D,rep,2);

% Permutation p-value
subplot(s,t,5)
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
%scatter(x,p(:,2), 500,'k.','jitter','on', 'jitterAmount', 0.3);
p=tN(:,k,l);
[f1,xi1]=ksdensity(p,'support',[-1,1]);
p=tN(:,n,n);
[f,xi]=ksdensity(p,'support',[-1,1]);
hold on
%plot(xi,f,'.-','LineWidth',2);
plot(xi,f,'.:',xi1,f1,'k.-','LineWidth',3);
set(gca,'FontSize',fontSize);
h=legend('Mcorr','MGC','Location','NorthEast');
set(h,'FontSize',fontSize);
plot(tA(end),0.1,'.',tA(k,l),0.1,'kx','MarkerSize',mkSize);

x1 = tA(end);
y1 = 0.01;
x2 = tA(k,l);
y2 = 0.01;
txt1 = strcat('Mcorr p = ',num2str((floor(1-pAll(end))*100)/100));
txt2 = strcat('MGC p < 1/',num2str(repp));
text(x1,y1,txt1)
text(x2,y2,txt2)


xlim([minp,maxp]);
%ylabel('Density Function','FontSize',fontSize);
title('Density Plot','FontSize',fontSize);
%ylim([-1 15]);
%set(gca,'YTick',[]); % Remove y axis ticks
hold off

% figure
% p=tN(:,k,l);
% %scatter(x,p(:,2), 500,'k.','jitter','on', 'jitterAmount', 0.3);
% [f1,xi1]=ksdensity(p,'support',[-1,1]);
% hold on
% plot(xi1,f1,'.-','LineWidth',2);
% plot(tA(k,l),0.01,'k.','MarkerSize',30);
% xlim([minp,maxp]);
% %ylim([-1 15]);
% set(gca,'FontSize',16);
% set(gca,'YTick',[]); % Remove y axis ticks
% hold off

%xlabel('False Positive Rate','FontSize',20);
%ylabel('Density Function','FontSize',20);
%title('Brain Activity vs. Fake Movie','FontSize',24);
1-pAll(end)
1-pAll(k,l)
power1(end)
power1(k,l)

% A2: heatmaps of the distance matrices
ax=subplot(s,t,2);
if samebar==1
    maxC=ceil(max(max([C,D]))*10)/10;
else
    maxC=ceil(max(max([C]))*10)/10;
end
imagesc(C');
colormap(ax,map3);
caxis([0,maxC]);
if colorb~=0
    h=colorbar('Ticks',[-1000,0,maxC]);
end
%xlabel('$$\tilde{A}_{ij}=|x_{i}-x_{j}|_{2}$$','FontSize',30,'interpreter','latex','position',[0,-1])
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
xlabel('$\|x_{i}-x_{j}\|$','interpreter','latex')
title('$\tilde{A}$','interpreter','latex');
pos = get(ax,'position');

%subplot(s,t,6)
ax=subplot(s,t,3);
if samebar~=1
    maxC=ceil(max(max([D]))*10)/10;
end
imagesc(D');
set(gca,'FontSize',16)
colormap(ax,map3);
if colorb~=0
    h=colorbar('Ticks',[-1000,0,maxC]);
end
caxis([0,maxC]);
% title('Distance Heatmap of Y','FontSize',16)
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
%xlabel('\tilde{B}_{ij}=\|y_{i}-y_{j}\|_{2}','FontSize',25);
xlabel('$\|y_{i}-y_{j}\|$','interpreter','latex')
title('$\tilde{B}$','interpreter','latex');

% A3: heatmaps of the doubly centered distance matrices
for option=1:1
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
    
    ax=subplot(s,t,6);
    if samebar==1
        minC=floor(min(min([A,B]))*10)/10;maxC=ceil(max(max([A,B]))*10)/10;
    else
        minC=floor(min(min([A]))*10)/10;maxC=ceil(max(max([A]))*10)/10;
    end
    minC=min(minC,-maxC);maxC=max(maxC,-minC);
    imagesc(A');
    colormap(ax,map2)
    caxis([minC,maxC]);
    if colorb~=0
        h=colorbar('Ticks',[-1000,minC,maxC]);
    end
    %title('Normalized Distances','FontSize',30)
    set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
    %xlabel('A=H\tilde{A}H','FontSize',25);
    title('$$A$$','FontSize',fontSize,'interpreter','latex');
    
    ax=subplot(s,t,7);
    if samebar==1
        minD=floor(min(min([A,B]))*10)/10;maxD=ceil(max(max([A,B]))*10)/10;
    else
        minD=floor(min(min([B]))*10)/10;maxD=ceil(max(max([B]))*10)/10;
    end
    minD=min(minD,-maxD);maxD=max(maxD,-minD);
    imagesc(B');
    set(gca,'FontSize',fontSize)
    colormap(ax,map2)
    caxis([minD,maxD]);
    if colorb~=0
        h=colorbar('Ticks',[-1000,minD,maxD]);
    end
    %title('Doubly-Centered Distance Heatmap of Y','FontSize',16)
    set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
    %xlabel('B=H\tilde{B}H','FontSize',25);
    title('$$B$$','interpreter','latex');
    
    % Local distance matrices
    %if n~=100 || noise~=1
    %   CorrIndTest(type,n,1,1,0, 1000,noise,0.05,[1,0,0,0]);
    %end
    %filename=strcat(pre1,'CorrIndTestType',num2str(type),'N',num2str(n),'Dim1.mat');
    %load(filename,'neighborhoods');
    %neighbor=neighborhoods(1,end);
    RC=(RC>k);
    RD=(RD>l);
    % ind=(dC.*dD==0);
    A(RC)=0;
    B(RD)=0;
    if cc==1
        A=A-mean(mean(A));B=B-mean(mean(B));
    end
    
    % A4: heatmaps of the local doubly centered distance matrices
    ax=subplot(s,t,9);
    imagesc(A');
    colormap(ax,map2)
    caxis([minC,maxC]);
    if colorb~=0
        h=colorbar('Ticks',[-1000,minC,maxC]);
    end
    %title('Local Normalized Distances','FontSize',30)
    %xlabel(strcat('k^{*}}=',num2str(k),'interpreter','latex'),'FontSize',25);
    title('$$A^{k^{*}}$$','interpreter','latex');
    set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
    
    ax=subplot(s,t,10);
    imagesc(B');
    set(gca,'FontSize',fontSize)
    colormap(ax,map2)
    caxis([minD,maxD]);
    if colorb~=0
        h=colorbar('Ticks',[-1000,minD,maxD]);
    end
    %title('Local Doubly-Centered Distance Heatmap of Y','FontSize',16)
    title('$$B^{l^{*}}$$','interpreter','latex');
    set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
    
    % A4: heatmaps of the distance covariance entries
    ax=subplot(s,t,4);
    MH=ceil(max(max(mcorrH(2:end,2:end)))*10)/10;
    mH=floor(min(min(mcorrH(2:end,2:end)))*10)/10;
    mH=min(mH,-MH);MH=max(MH,-mH);
    imagesc(mcorrH');
    colormap(ax,map2)
%     if colorb~=0
        h=colorbar('Ticks',[mH,0,MH]);
%     end
    caxis([mH,MH]);
    %ylabel('$$A^{k^{*}}.*B^{l^{*}}$$','FontSize',25,'Rotation',0,'position',[-5,50],'interpreter','latex');
    set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
    %title('Element-wise Product of','FontSize',30)
    title('$$A.*B$$','FontSize',fontSize,'interpreter','latex');
    pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
    %title('Local Distance Covariance of (X, Y)','FontSize',16)
    
    % A4: heatmaps of the distance covariance entries
    ax=subplot(s,t,8);
    mcorrH=A.*B;
    %     mcorrH(RC)=0;
    %     mcorrH(RD)=0;
    %     %MH=ceil(max(max(mcorrH(2:end,2:end)))*10)/10;
    %     MH=ceil(max(max(mcorrH(2:end,2:end)))*10)/10;
    %     mMH=floor(min(min(mcorrH(2:end,2:end)))*10)/10;
    imagesc(mcorrH');
    colormap(ax,map2)
    if colorb~=0
        h=colorbar('Ticks',[-1000,mH,MH]);
    end
    caxis([mH,MH]);
    %ylabel('$$A^{k^{*}}.*B^{l^{*}}$$','FontSize',25,'Rotation',0,'position',[-5,50],'interpreter','latex');
    set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove y axis ticks
    title('$$A^{k^{*}}.*B^{l^{*}}$$','interpreter','latex');
    %title('Element-wise Product of','FontSize',30)
    %title('Local Distance Covariance of (X, Y)','FontSize',16)
    
    ax=subplot(s,t,12);
    %filename=strcat(pre1,'CorrIndTestType',num2str(type),'N100Dim1.mat');
    %load(filename)
    kmin=2;
    hold on
    ph=power1(kmin:n,kmin:n)';
%     ph(ph>0.98)=0.98;
    %     nei=find(ph>0.99);
    %     [k,l]=ind2sub(size(ph),nei);
    imagesc(ph);
    set(gca,'FontSize',15)
    set(gca,'YDir','normal')
    cmap=map4;
%     cmap=zeros(50,3);
%     cmap(1:45,:)=map4(1:45,:);
%     cmap(46:49,:)=map4(80:83,:);
%     cmap(50,:)=map4(128,:);
    colormap(ax,cmap)
    caxis([0 1])
    xlabel('k','FontSize',fontSize)
    ylabel('l','FontSize',fontSize,'Rotation',0,'position',[-7,20]);
    xlim([1 n-1]);
    ylim([1 n-1]);
    h=colorbar('Ticks',[0,0.5,1]);
    set(h,'FontSize',fontSize);
        pos2 = get(ax,'position');
        pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
    
%     imagesc(k,l,1);
    
    hold off
    %     %set(h,'FontSize',16,'location','southoutside');
    %set(gca,'XTick',[]); % Remove x axis ticks
    %set(gca,'YTick',[]); % Remove y axis ticks
    title('Multiscale Power Map','FontSize',fontSize);
    colorbar
end

%
F.fname=strcat(pre2, 'A2');
F.wh=[8 4]*2;
print_fig(gcf,F)
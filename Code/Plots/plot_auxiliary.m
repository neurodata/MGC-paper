function []=plot_auxiliary(type,pre1)%,pre2)

% type=6;n=100;dim=1;noise=1;pre1='../../Data/';
% CorrSimPlotsA(type,n,dim,noise,pre1);
% Used to generate figure A in the files

if nargin<1
    type=11;
end
if nargin<2
    pre1='../../Data/Results/'; % The folder to locate data
end
% if nargin<2
%     pre2='../../Draft/Figures/FigReal'; % The folder to save figures
% end
option=2;
n=50;
dim=1;
noise=0;
cc=1;
samebar=1;
colorb=1;
rep=1000;

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
set(groot,'defaultAxesColorOrder',map1);

% Generate data
[x, y]=CorrSampleGenerator(type,n,dim,1, noise);
if noise~=0
[x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
end
% A1: scatter plot
figure
if noise==0
plot(x,y,'.','MarkerSize',30);
else
plot(x,y,'.',x1,y1,'k.','MarkerSize',30);
end
xlabel('X','FontSize',25)
ylabel('Y','FontSize',25,'Rotation',0,'position',[-1.3,-0.1]);
if samebar==1
    xlim([-1.2,1.2]);
    ylim([-1.2,1.2]);
end
title('Scatter Plot of (X,Y)','FontSize',30)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
% Get distance matrix
[x,ind]=sort(x,'ascend');
y=y(ind);
C=squareform(pdist(x));
D=squareform(pdist(y));

% optimal scale
tA=zeros(n,n,rep);
tN=zeros(n,n,rep);
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
% figure
% ksdensity(reshape(tA(n,n,:),1,rep));
% figure
% ksdensity(reshape(tN(n,n,:),1,rep));
power1(1,:)=0;power1(:,1)=0; % Set the powers of all local tests at rank 0 to 0

neighbor=maxNeighbors(power1,0,tN,tA);
% [~,pAll]=MGCPermutationTest(C,D,rep,2);

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
figure
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
plot(xi,f,'.:',xi1,f1,'k.-','LineWidth',2);
set(gca,'FontSize',24);
h=legend('Mcorr','MGC','Location','NorthEast');
set(h,'FontSize',30);
plot(tA(end),0.01,'.',tA(k,l),0.01,'k.','MarkerSize',30);
xlim([minp,maxp]);
ylabel('Density Function','FontSize',24);
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
% 1-pAll2(end)
% 1-pAll2(k,l)

% A2: heatmaps of the distance matrices
if samebar==1
    maxC=ceil(max(max([C,D]))*10)/10;
else
    maxC=ceil(max(max([C]))*10)/10;
end
figure
imagesc(C');
colormap(map3)
caxis([0,maxC]);
if colorb~=0
    h=colorbar('Ticks',[-1000,0,maxC]);
    set(h,'FontSize',16);
end
title('Pairwise Distances','FontSize',30)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
%xlabel('\tilde{A}_{ij}=\|x_{i}-x_{j}\|_{2}','FontSize',25);
ylabel('$$\tilde{A}$$','FontSize',30,'Rotation',0,'position',[-5,50],'interpreter','latex');

%subplot(s,t,6)
figure
if samebar~=1
    maxC=ceil(max(max([D]))*10)/10;
end
imagesc(D');
set(gca,'FontSize',16)
colormap(map3)
if colorb~=0
    h=colorbar('Ticks',[-1000,0,maxC]);
    set(h,'FontSize',16);
end
caxis([0,maxC]);
% title('Distance Heatmap of Y','FontSize',16)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
%xlabel('\tilde{B}_{ij}=\|y_{i}-y_{j}\|_{2}','FontSize',25);
ylabel('$$\tilde{B}$$','FontSize',30,'Rotation',0,'position',[-5,50],'interpreter','latex');

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
    
    if samebar==1
        minC=floor(min(min([A,B]))*10)/10;maxC=ceil(max(max([A,B]))*10)/10;
    else
        minC=floor(min(min([A]))*10)/10;maxC=ceil(max(max([A]))*10)/10;
    end
    figure
    imagesc(A');
    colormap(map2)
    caxis([minC,maxC]);
    if colorb~=0
        h=colorbar('Ticks',[-1000,minC,maxC]);
        set(h,'FontSize',16);
    end
    %title('Normalized Distances','FontSize',30)
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    %xlabel('A=H\tilde{A}H','FontSize',25);
    ylabel('$$A$$','FontSize',30,'Rotation',0,'position',[-5,50],'interpreter','latex');
    
    figure
    if samebar==1
        minD=floor(min(min([A,B]))*10)/10;maxD=ceil(max(max([A,B]))*10)/10;
    else
        minD=floor(min(min([B]))*10)/10;maxD=ceil(max(max([B]))*10)/10;
    end
    imagesc(B');
    set(gca,'FontSize',16)
    colormap(map2)
    caxis([minD,maxD]);
    if colorb~=0
        h=colorbar('Ticks',[-1000,minD,maxD]);
        set(h,'FontSize',16);
    end
    %title('Doubly-Centered Distance Heatmap of Y','FontSize',16)
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    %xlabel('B=H\tilde{B}H','FontSize',25);
    ylabel('$$B$$','FontSize',30,'Rotation',0,'position',[-5,50],'interpreter','latex');
    
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
    figure
    imagesc(A');
    colormap(map2)
    caxis([minC,maxC]);
    if colorb~=0
        h=colorbar('Ticks',[-1000,minC,maxC]);
        set(h,'FontSize',16);
    end
    %title('Local Normalized Distances','FontSize',30)
    %xlabel(strcat('k^{*}}=',num2str(k),'interpreter','latex'),'FontSize',25);
    title('$$A^{k^{*}}$$','FontSize',30,'interpreter','latex');
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    
    figure
    imagesc(B');
    set(gca,'FontSize',16)
    colormap(map2)
    caxis([minD,maxD]);
    if colorb~=0
        h=colorbar('Ticks',[-1000,minD,maxD]);
        set(h,'FontSize',16);
    end
    %title('Local Doubly-Centered Distance Heatmap of Y','FontSize',16)
    title('$$B^{l^{*}}$$','FontSize',30,'interpreter','latex');
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    
    % A4: heatmaps of the distance covariance entries
    figure
    MH=ceil(max(max(mcorrH(2:end,2:end)))*10)/10;
    mMH=floor(min(min(mcorrH(2:end,2:end)))*10)/10;
    imagesc(mcorrH');
    set(gca,'FontSize',16)
    colormap(map2)
    if colorb~=0
        h=colorbar('Ticks',[-1000,mMH,MH]);
        set(h,'FontSize',16);
    end
    caxis([mMH,MH]);
    %ylabel('$$A^{k^{*}}.*B^{l^{*}}$$','FontSize',25,'Rotation',0,'position',[-5,50],'interpreter','latex');
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    %title('Element-wise Product of','FontSize',30)
    title('$$A.*B$$','FontSize',30,'interpreter','latex');
    %title('Local Distance Covariance of (X, Y)','FontSize',16)
    
    % A4: heatmaps of the distance covariance entries
    figure
    mcorrH=A.*B;
    %     mcorrH(RC)=0;
    %     mcorrH(RD)=0;
    %     %MH=ceil(max(max(mcorrH(2:end,2:end)))*10)/10;
    %     MH=ceil(max(max(mcorrH(2:end,2:end)))*10)/10;
    %     mMH=floor(min(min(mcorrH(2:end,2:end)))*10)/10;
    imagesc(mcorrH');
    set(gca,'FontSize',16)
    colormap(map2)
    if colorb~=0
        h=colorbar('Ticks',[-1000,mMH,MH]);
        set(h,'FontSize',16);
    end
    caxis([mMH,MH]);
     if colorb~=0
        h=colorbar('Ticks',[-1000,mMH,MH]);
        set(h,'FontSize',16);
    end
    %ylabel('$$A^{k^{*}}.*B^{l^{*}}$$','FontSize',25,'Rotation',0,'position',[-5,50],'interpreter','latex');
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    title('$$A^{k^{*}}.*B^{l^{*}}$$','FontSize',30,'interpreter','latex');
    %title('Element-wise Product of','FontSize',30)
    %title('Local Distance Covariance of (X, Y)','FontSize',16)
    
    figure
    %filename=strcat(pre1,'CorrIndTestType',num2str(type),'N100Dim1.mat');
    %load(filename)
    kmin=2;
    hold on
    ph=power1(kmin:n,kmin:n)';
    ph(ph>0.98)=0.98;
    %     nei=find(ph>0.99);
    %     [k,l]=ind2sub(size(ph),nei);
    imagesc(ph);
    set(gca,'FontSize',16)
    set(gca,'YDir','normal')
    map2 = brewermap(128,'GnBu'); % brewmap
    cmap=zeros(50,3);
    cmap(1:48,:)=map2(1:48,:);
    cmap(49,:)=map2(80,:);
    cmap(50,:)=map2(128,:);
    colormap(cmap)
    caxis([0.5 1])
    xlabel('k','FontSize',25)
    ylabel('l','FontSize',25)%,'Rotation',0,'position',[-5,50]);
    xlim([1 n-1]);
    ylim([1 n-1]);
    h=colorbar('Ticks',[-1000,0.5,0.98]);
    set(h,'FontSize',16);
    
    imagesc(k,l,1);
    
    hold off
    %     %set(h,'FontSize',16,'location','southoutside');
    %set(gca,'XTick',[]); % Remove x axis ticks
    %set(gca,'YTick',[]); % Remove y axis ticks
    title('Testing Powers','FontSize',30);
end
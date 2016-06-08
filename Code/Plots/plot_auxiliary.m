function []=CorrSimPlotsA(type,n,dim,noise,pre1,pre2)

% type=6;n=100;dim=1;noise=1;pre1='../../Data/'; 
% CorrSimPlotsA(type,n,dim,noise,pre1);
% Used to generate figure A in the files

if nargin<1
    type=6; 
end
if nargin<2
    n=100; 
end
if nargin<3
    dim=1; 
end
if nargin<4
    noise=1; 
end
if nargin<5
    pre1='../../Data/Results/'; % The folder to locate data
end
if nargin<6
    pre2='../../../../Draft/Figures/FigReal'; % The folder to save figures
end

cmap=zeros(3,3);
gr = [0.5,0.5,0.5];
ma = [1,0,1];
% cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = gr;
% cmap(3,:) = cy;
map1=cmap;
map2 = brewermap(128,'BuGn'); % brewmap
set(groot,'defaultAxesColorOrder',map1);

% Generate data
[x, y]=CorrSampleGenerator(type,n,dim,1, noise);
[x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
% A1: scatter plot
figure
plot(x,y,'.',x1,y1,'k.','MarkerSize',30);
xlabel('X','FontSize',25)
ylabel('Y','FontSize',25,'Rotation',0,'position',[-1.1,0.25]);
title('Scatter Plot of (X,Y)','FontSize',30)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
        
% Get distance matrix
[x,ind]=sort(x,'ascend');
y=y(ind);
C=squareform(pdist(x));
D=squareform(pdist(y));
dC=disToRanks(C);
dD=disToRanks(D);
H=eye(n)-(1/n)*ones(n,n);
A=H*C;
B=D*H;
% % For mcorr, further adjust the double centered matrices to remove high-dimensional bias
% A=A-C/n;
% B=B-D/n;
% for j=1:n
%     A(j,j)=0;
%     B(j,j)=0;
% end

%global stat
mcorrH=A.*B;

% A2: heatmaps of the distance matrices
maxC=max(max(C));
maxD=max(max(D));
maxC=ceil(max(maxC,maxD));
figure
imagesc(C');
colormap(map2)
caxis([0,maxC]);
title('Pairwise Distances','FontSize',30)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
xlabel('\tilde{A}_{ij}=\|x_{i}-x_{j}\|_{2}','FontSize',25);
ylabel('X','FontSize',25,'Rotation',0,'position',[-5,50]);
%subplot(s,t,6)
figure
imagesc(D');
set(gca,'FontSize',16)
colormap(map2)
h=colorbar('Ticks',[-1,0,maxC]);
set(h,'FontSize',16,'location','southoutside');
caxis([0,maxC]);
% title('Distance Heatmap of Y','FontSize',16)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
xlabel('\tilde{B}_{ij}=\|y_{i}-y_{j}\|_{2}','FontSize',25);
ylabel('Y','FontSize',25,'Rotation',0,'position',[-5,50]);

% A3: heatmaps of the doubly centered distance matrices
minC=min(min([A,B]));maxC=max(max([A,B]));
minD=floor(minC);maxD=ceil(maxC);
C=A;
D=B;
figure
imagesc(C');
colormap(map2)
caxis([minD,maxD]);
title('Normalized Distances','FontSize',30)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
xlabel('A=H\tilde{A}H','FontSize',25);
ylabel('X','FontSize',25,'Rotation',0,'position',[-5,50]);
figure
imagesc(D');
set(gca,'FontSize',16)
colormap(map2)
caxis([minD,maxD]);
h=colorbar('Ticks',[-1000,minD,maxD]);
set(h,'FontSize',16,'location','southoutside');
%title('Doubly-Centered Distance Heatmap of Y','FontSize',16)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
xlabel('B=H\tilde{B}H','FontSize',25);
ylabel('Y','FontSize',25,'Rotation',0,'position',[-5,50]);

% Local distance matrices
if n~=100 || noise~=1
    CorrIndTest(type,n,1,1,0, 1000,noise,0.05,[1,0,0,0]);
end
filename=strcat(pre1,'CorrIndTestType',num2str(type),'N',num2str(n),'Dim1.mat');
load(filename,'neighborhoods');
neighbor=neighborhoods(1,end);
l=ceil(neighbor/n);
k=neighbor-n*(l-1);
dC=(dC>k);
dD=(dD>l);
% ind=(dC.*dD==0);
A(dC)=0;
B(dD)=0;

% A4: heatmaps of the local doubly centered distance matrices
C(dC)=minD;
D(dD)=minD;
figure
imagesc(C');
colormap(map2)
caxis([minD,maxD]);
title('Local Normalized Distances','FontSize',30)
xlabel('A^{k^{*}}=','FontSize',25);
ylabel('X','FontSize',25,'Rotation',0,'position',[-5,50]);
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
figure
imagesc(D');
set(gca,'FontSize',16)
colormap(map2)
caxis([minD,maxD]);
h=colorbar('Ticks',[-1000,minD,maxD]);
set(h,'FontSize',16,'location','southoutside');
%title('Local Doubly-Centered Distance Heatmap of Y','FontSize',16)
ylabel('Y','FontSize',25,'Rotation',0,'position',[-5,50]);
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks

% A4: heatmaps of the distance covariance entries
figure
mcorrH(dC)=0;
mcorrH(dD)=0;
MH=ceil(max(max(mcorrH(2:end,2:end))));
imagesc(mcorrH');
set(gca,'FontSize',16)
colormap(map2)
%h=colorbar;
%set(h,'FontSize',16);
caxis([0,MH]);
h=colorbar('Ticks',[-1000,0,MH]);
set(h,'FontSize',16,'location','southoutside');
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
title('Element-wise Product of','FontSize',30)
%title('Local Distance Covariance of (X, Y)','FontSize',16)
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
    pre1='../../Data/'; % The folder to locate data
end
if nargin<6
    pre2='../../Figures/Fig'; % The folder to save figures
end

% Generate data
[x, y]=CorrSampleGenerator(type,n,dim,1, noise);
[x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
% A1: scatter plot
figure
plot(x,y,'b.',x1,y1,'r.');
xlabel('X','FontSize',16)
ylabel('Y','FontSize',16);
title('Scatter Plot of (X,Y)','FontSize',16)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
        
% Get distance matrix
[x,ind]=sort(x,'ascend');
y=y(ind);
C=squareform(pdist(x));
D=squareform(pdist(y));
dC=disToRanks(C);
dD=disToRanks(D);
map2 = brewermap(128,'GnBu'); % brewmap
H=eye(n)-(1/n)*ones(n,n);
A=H*C*H;
B=H*D*H;
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
maxC=max(maxC,maxD);
figure
imagesc(C');
colormap(map2)
caxis([0,maxC]);
title('Distance Heatmap of X','FontSize',16)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
%subplot(s,t,6)
figure
imagesc(D');
set(gca,'FontSize',16)
colormap(map2)
h=colorbar;
set(h,'FontSize',16);
caxis([0,maxC]);
title('Distance Heatmap of Y','FontSize',16)

% A3: heatmaps of the doubly centered distance matrices
minC=min(min([A,B]));maxC=max(max([A,B]));
minD=minC;maxD=maxC;
C=A;
D=B;
figure
imagesc(C');
colormap(map2)
caxis([minC,maxC]);
title('Doubly-Centered Distance Heatmap of X','FontSize',16)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
figure
imagesc(D');
set(gca,'FontSize',16)
colormap(map2)
h=colorbar;
set(h,'FontSize',16);
caxis([minD,maxD]);
title('Doubly-Centered Distance Heatmap of Y','FontSize',16)

% Local distance matrices
if n~=100 || noise~=1
    CorrIndTest(type,n,1,1,0, 1000,noise,0.05,[1,0,0,0]);
end
filename=strcat(pre1,'CorrIndTestType',num2str(type),'N',num2str(n),'Dim1.mat');
load(filename,'neighborhoods');
neighbor=neighborhoods(1,end);
l=ceil(neighbor/n);
k=neighbor-n*(l-1);
dC=(dC<k);
dD=(dD<l);
ind=(dC.*dD==0);
A(ind)=0;
B(ind)=0;

% A4: heatmaps of the local doubly centered distance matrices
C(ind)=minC;
D(ind)=minD;
figure
imagesc(C');
colormap(map2)
caxis([minC,maxC]);
title('Local Doubly-Centered Distance Heatmap of X','FontSize',16)
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
figure
imagesc(D');
set(gca,'FontSize',16)
colormap(map2)
caxis([minD,maxD]);
title('Local Doubly-Centered Distance Heatmap of Y','FontSize',16)

% A4: heatmaps of the distance covariance entries
figure
mcorrH(ind)=0;
imagesc(mcorrH');
set(gca,'FontSize',16)
colormap(map2)
h=colorbar;
set(h,'FontSize',16);
caxis([0,1]);
title('Local Distance Covariance of (X, Y)','FontSize',16)
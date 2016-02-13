function []=CorrSimPlotsA(type,n,dim,noise,pre1,pre2)

% type=6;n=20;dim=1;noise=0;pre1='../../Data/'; 
% CorrSimPlotsA(type,n,dim,noise,pre1);

if nargin<1
    type=6; % Usually 20, but can be changed in case of new simulations
end
if nargin<2
    n=20; % Usually 20, but can be changed in case of new simulations
end
if nargin<3
    dim=1; % Usually 20, but can be changed in case of new simulations
end
if nargin<4
    noise=0; % Usually 20, but can be changed in case of new simulations
end
if nargin<5
    pre1='../../Data/'; % The folder to locate data
end
if nargin<6
    pre2='../../Figures/Fig'; % The folder to save figures
end

%%%performance profile
figNumber='A';
%figure('units','normalized','position',[0 0 1 1])
%set(groot,'defaultAxesColorOrder',map1);
s=2;
t=2;
tryy=2;

% A1
[x, y]=CorrSampleGenerator(type,n,dim,1, noise);
%[x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
% figure
% plot(x,y,'b.',x1,y1,'r.');
        
% A2
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
% For mcorr, further adjust the double centered matrices to remove high-dimensional bias
A=A-C/n;
B=B-D/n;
for j=1:n
    A(j,j)=0;
    B(j,j)=0;
end

%global stat
mcorr=sum(sum(A.*B))/norm(A,'fro')/norm(B,'fro');

if tryy==1
   minC=min(min([A,B]));maxC=max(max([A,B]));
   minD=minC;maxD=maxC;
   C=A;
   D=B;
else
   maxC=max(max(C));
   maxD=max(max(D));
   minC=0;
   minD=0;
end
subplot(s,t,1)
imagesc(C');
colormap(map2)
colorbar
caxis([minC,maxC]);
title('Global Distance Heatmap of X')
subplot(s,t,2)
imagesc(D');
colormap(map2)
colorbar
caxis([minD,maxD]);
title('Global Distance Heatmap of Y')

% A3
if n~=100 || noise~=1
    CorrIndTest(type,n,1,1,0, 1000,noise,0.05,[1,0,0,0]);
end
filename=strcat(pre1,'CorrIndTestType',num2str(type),'N',num2str(n),'Dim1.mat');
load(filename,'power1All');
neighbor=verifyNeighbors(1-power1All(:,:,end));
l=ceil(neighbor/n);
k=neighbor-n*(l-1);
dC=(dC<k);
dD=(dD<l);
ind=(dC.*dD==0);
A(ind)=0;
B(ind)=0;
% figure
% gplot(1-ind,[x,y],':or');
% h=findobj('type','line');
% set(h,'MarkerSize',10)

% A4
C(ind)=0;
D(ind)=0;
subplot(s,t,3)
imagesc(C');
colormap(map2)
colorbar
caxis([minC,maxC]);
title('Local Distance Heatmap of X')
subplot(s,t,4)
imagesc(D');
colormap(map2)
colorbar
caxis([minD,maxD]);
title('Local Distance Heatmap of Y')

%local stat
lmcorr=sum(sum(A.*B))/norm(A,'fro')/norm(B,'fro');

%F.fname=strcat(pre2, figNumber);
%F.wh=[3 2.5]*2;
%print_fig(gcf,F)
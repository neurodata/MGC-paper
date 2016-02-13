function [ppower1,ppower2]=CorrSimPermScale(n,dim,rep1,rep2,ratio,noise,alpha)
% Author: Cencheng Shen
% n=100;dim=1;rep1=1;rep2=100;ratio=1;noise=1;
% [p1,p2]=CorrSimPermScale(n,dim,rep1,rep2,ratio,noise);
% Used to plot figure 0 in the files
if nargin<1
    n=100;
end
if nargin<2
    dim=1;
end
if nargin<3
    rep1=10;
end
if nargin<4
    rep2=1000; % The folder to save figures
end
if nargin<5
    ratio=10;
end
if nargin<6
    noise=1;
end
if nargin<7
    alpha=0.05; % Default type 1 error level
end

option=[1,0,0];
total=20;
%figure('units','normalized','position',[0 0 1 1])
if dim>1
    noise=0;
end

pre1='../../Data/'; 
pre2='../../Figures/Fig'; % The folder to save figures
ppower1=zeros(n,n,total);
ppower2=zeros(n,n,total);
for type=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
    load(filename,'power1All');
    ppower2(:,:,type)=power1All(:,:,end);
    for r1=1:rep1
        [x, y]=CorrSampleGenerator(type,n,dim,1, noise);
        [p1]=CorrPermDistTest(squareform(pdist(x)),squareform(pdist(y)),rep2,0,'Ind',ratio,zeros(3,1),alpha,option);
        ppower1(:,:,type)=ppower1(:,:,type)+p1/rep1;
    end
end

map2 = brewermap(128,'GnBu'); % brewmap
figNumber='PermSim';
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    K=n;kmin=1;
    ph=ppower1(kmin:K,kmin:K,j)';
    imagesc(ph);
    set(gca,'YDir','normal')
    colormap(map2)
    %caxis([0 thres])
    colorbar();
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    title(titlechar);
end
xlabel('Neighborhood Choice of X','position',[-200 -20],'FontSize',20);
ylabel('Neighborhood Choice of Y','position',[-535 300],'FontSize',20);
colorbar
tstring=' of mcorr ';
h=suptitle(strcat('Testing Powers of All Local Tests',tstring, ' for Dimension 1'));
set(h,'FontSize',20,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 4]*2;
print_fig(gcf,F)

figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    K=n;kmin=1;
    ph=ppower2(kmin:K,kmin:K,j)';
    imagesc(ph);
    set(gca,'YDir','normal')
    colormap(map2)
    %caxis([0 thres])
    colorbar();
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    title(titlechar);
end
xlabel('Neighborhood Choice of X','position',[-200 -20],'FontSize',20);
ylabel('Neighborhood Choice of Y','position',[-535 300],'FontSize',20);
colorbar
tstring=' of mcorr ';
h=suptitle(strcat('Testing Powers of All Local Tests',tstring, ' for Dimension 1'));
set(h,'FontSize',20,'FontWeight','normal');

filename=strcat(pre1,'CorrSimPermScale',num2str(type),'N',num2str(n),'Dim',num2str(dim));
save(filename,'ppower1','ppower2','n','rep1','rep2','dim','noise','alpha','option');
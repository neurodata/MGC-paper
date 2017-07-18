%% simulate data
type=8; % simulation type, can be set to any of 1-20
d=1; % dimension, recommend to be 1 for 2D visualization
noise=0; % noise level, can be set to any value no smaller than 0
n=100; % sample size, 

[x,y]=CorrSampleGenerator(type,n,d,1, noise);

%% run mgc
rep=1000;
option='mcor';

A=squareform(pdist(x));
B=squareform(pdist(y));

[pMGC,statMGC,pLocalCorr,localCorr,optimalInd]=MGCPermutationTest(A,B,rep,option);
display(strcat('The p-value and test statistic of MGC are',{' '},num2str(pMGC),{' '},'and',{' '},num2str(statMGC)));
display(strcat('For p-value, the smaller the more significant'));

%% figure output

figure
fs=15;
plot(x,y,'.','markerSize',15);
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
title('Visualization of the relationship','FontSize',fs)

figure
map2 = brewermap(128,'BuPu'); % brewmap
cmap2=flipud(map2);
cticks=[0.001, 0.01, 0.1, 0.5];
lw=1.5;
glob= [0.5,0.5,0.5];
set(groot,'defaultAxesColorOrder',cmap2);
hold on
imagesc(log((pLocalCorr+0.0001)'));
set(gca,'YDir','normal')
colormap(cmap2)
h=colorbar('Ticks',log(cticks),'TickLabels',cticks,'location','eastoutside','FontSize',fs-5);
set(gca,'FontSize',fs);
%     if i==3
xlabel('# Neighbors in X','FontSize',fs)
ylabel('# Neighbors in Y','FontSize',fs)
title('Significance Map and Optimal Scales','FontSize',fs)


%[~,indP]=MGCScaleVerify(p2All',rep);
indP=optimalInd;
[m,n]=size(pLocalCorr);
[J,I]=ind2sub([m,n],indP);
Ymin=min(I);
Ymax=max(I);
Xmin=min(J);
Xmax=max(J);
%
if Xmin==Xmax && Ymin==Ymax
    plot(Xmin,Ymin,'go','markerSize',6,'linewidth',3);
else
    plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
    plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
    plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
    plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
end
tmp=zeros(m,n);
tmp(J,I)=1;
tmp(localCorr<statMGC)=0;
[k,l]=ind2sub([m,n],find(tmp==1,1,'last'));
plot(k,l,'go','markerSize',6,'linewidth',3);
plot(m,n,'.','markerSize',18,'MarkerFaceColor',glob,'Color',glob)
xlim([1,m]);
ylim([1,n]);

function []=plot_realData
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
pre2=strcat(rootDir,'Figures/FigReal');% The folder to save figures

%% figure stuff

cmap=zeros(4,3);
gr =[0,1,0];
ma = [1,0,1];
map3 = brewermap(128,'PiYG'); % brewmap
lgr=map3(100,:);
dgr=map3(128,:);
gr=lgr;
glob= [0.5,0.5,0.5];
% cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = gr;
cmap(3,:) = gr;
cmap(4,:) = gr;
% cmap(3,:) = cy;
map1=cmap;
map2 = brewermap(128,'BuPu'); % brewmap
set(groot,'defaultAxesColorOrder',map1);


fnames={'CorrPermDistTestTypeBrainCxP.mat'; ...
    %'CorrPermDistTestTypeBrainLMRxY.mat'; ...
    'CorrPermDistTestTypeMigrainxCCI.mat'};
xlabs={ '# Activity Neighbors'; ...
       % '# Shape Neighbors';...
        '# Graph Neighbors'};
ylabs={ '# Personality Neighbors'; ...
       % '# Disease Neighbors';...
        '# Creativity Neighbors'};
% tits= {'A. Brain Activity vs. Personality'; ...
%     'B. Brain Shape vs. Disorder';...
%     'C. Brain Graph vs. Creativity'};
tits={'i';'ii';'i';'ii';'i';'ii'};
s=2;t=3;

height=0.3;
vspace=0.15;

width=0.27;
left1=0.08;
left2=0.35;
left3=0.65;
%left=width+left1;
bottom=nan(1,4);
bottom1=0.6;
bottom2=0.1;
for i=1:4
    bottom(i)=bottom1+(i-1)*(height+vspace);
end
bottom(4)=bottom(4)+0.15;
bottom(3)=bottom(3)+0.15;
bottom(2)=bottom(2)+0.08;
% bottom(1)=0.05;

F.pos1 =[left1, bottom1, width, height];
F.pos2=[left1, bottom2, width, height];
F.pos3=[left2, bottom1, width, height];
F.pos4=[left2, bottom2, width, height];
F.pos5=[left3, bottom1, width, height];
F.pos6=[left3, bottom2, width, height];


%% loop maps
figure(1), clf, hold all
filename=strcat(pre1,fnames{2});
load(filename);
cmap2=flipud(map2);
fs=9;
cticks=[0.001, 0.01, 0.1, 0.5];
lw=1.5;

for i=1:2
    filename=strcat(pre1,fnames{i});
    load(filename);
    [m,n]=size(pMLocal);
    
    if i==1
        subplot('Position',F.pos1)
    else
        subplot('Position',F.pos2)
    end
    hold on
    imagesc(log(pMLocal'));
    set(gca,'YDir','normal')
    colormap(cmap2)
    set(gca,'FontSize',fs);
    %     if i==3
    xlabel(xlabs{i},'FontSize',fs, ...
       'Units', 'normalized','Position', [-0.01, -0.16], 'HorizontalAlignment', 'left')
    ylabel(ylabs{i},'FontSize',fs, ...
        'Units', 'normalized','Position', [-0.21 0], 'HorizontalAlignment', 'left')
    title(strcat(tits{i},'.',{' '},'Significance Map'),'FontSize',fs, ...
        'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')
    
    
    %[~,indP]=MGCScaleVerify(p2All',rep);
    indP=optimalInd;
    [J,I]=ind2sub([m,n],indP);
    Ymin=min(I);
    Ymax=max(I);
    Xmin=min(J);
    Xmax=max(J);
    %
    if Xmin==Xmax && Ymin==Ymax
         plot(Xmin,Ymin,'g.','markerSize',20);
    else
        plot([Xmin,Xmin],[Ymin,Ymax],'g','linewidth',lw)
        plot([Xmax,Xmax],[Ymin,Ymax],'g','linewidth',lw)
        plot([Xmin,Xmax],[Ymin,Ymin],'g','linewidth',lw)
        plot([Xmin,Xmax],[Ymax,Ymax],'g','linewidth',lw)
    end
    tmp=zeros(m,n);
    tmp(J,I)=1;
    tmp(testMLocal<testMGC)=0;
    [k,l]=ind2sub([m,n],find(tmp==1,1,'last'));
    plot(k,l,'g.','markerSize',20);
    plot(m,n,'.','markerSize',20,'MarkerFaceColor',glob,'Color',glob)
    xticks=[5,round(m/2)-1,m-1];
    if i==1,  xticks(1)=3; end
    %  set(gca,'XTick',xticks,'XTickLabel',[2,round(m/2),m]); % Remove x axis ticks
%    set(gca,'YTick',[3,round(n/2)-1,n-1],'YTickLabel',[2,round(n/2),n]); % Remove x axis ticks
    xlim([2,m]);
    ylim([2,n]);
    axis('square');
    hold off
    if i==1
        h=colorbar('Ticks',log(cticks),'TickLabels',cticks,'location','eastoutside','FontSize',fs-3);
        title(h,'p-value')
    end
end

% plot last figure
load(strcat(pre1,'CorrBrainNoiseSummary.mat'));
% cmap=zeros(3,3);
% ma = [1,0,1];
% cmap(1,:) = ma;
% cmap(2,:) = ma;
% map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

subplot('Position',F.pos6)
%scatter(x,p(:,2), 500,'k.','jitter','on', 'jitterAmount', 0.3);
pv=p(:,1);
[f,xi]=ksdensity(pv);
hold on
plot(xi,f,'.-','LineWidth',lw);
pv=sort(pv,'ascend');
ord=0.01*ones(length(pv),1);
for i=2:length(pv);
    if pv(i)-pv(i-1)<0.001
        ord(i)=ord(i-1)+0.4;
    end
end
plot(pv,ord,'.','MarkerSize',8);
% <<<<<<< HEAD
xlim([0,0.15]);
ylim([-1 max(f)+1]);
% =======
% xlim([-0.05,0.15]);
% ylim([-1 15]);
% >>>>>>> 7132b2753b089edc2ca608140c119b901e31b17a
set(gca,'FontSize',fs);
set(gca,'YTick',[]); % Remove y axis ticks
axis('square');
hold off
xlabel('False Positive Rate','FontSize',fs,...
    'Units', 'normalized','Position', [-0.01, -0.16], 'HorizontalAlignment', 'left')
ylabel('Density','FontSize',fs, ...
    'Units', 'normalized','Position', [-0.15 0], 'HorizontalAlignment', 'left')
% title('D. Brain Activity vs. Fake Movie','FontSize',fs, ...
%     'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')
title(strcat(tits{6},'.',{' '},'Neuroimaging FPR'),'FontSize',fs, ...
    'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left')

% F.fname=strcat(pre2, 'CORR');
% F.wh=[3 2.5]*2;
% print_fig(gcf,F)
%colorbar()

load(strcat(rootDir,'Data/Preprocessed/BrainHippoShape.mat'));
subplot('Position',F.pos4)
% [A,B,RX,RY]=MGCDistTransform(LMRS,LMLS);
% k=60;
% LMRS(RX>k)=0;
per1=(Label==1);
per2=(Label==2);
per3=(Label==3);
t1=LMRS(per2,per2);
t2=LMRS(per2,per1);
t3=LMRS(per2,per3);
t1=t1(t1>0);t2=t2(t2>0);t3=t3(t3>0);
t2=[t2;t3];

pv=t1;
[f,xi]=ksdensity(pv);
hold on
h1=plot(xi,f,'b.-','LineWidth',2);

pv=t2;
[f,xi]=ksdensity(pv);
h2=plot(xi,f,'r.-','LineWidth',2);
xlim([0,max(xi)]);
hold off
ylabel('Density','FontSize',fs, ...
    'Units', 'normalized','Position', [-0.15 0], 'HorizontalAlignment', 'left')
xlabel('Distance','FontSize',fs,...
    'Units', 'normalized','Position', [-0.01, -0.16], 'HorizontalAlignment', 'left')
set(gca,'FontSize',fs);
title(strcat(tits{4},'.',{' '},'Distance Density Plot'),'FontSize',fs, ...
    'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left');
set(gca,'YTick',[],'YTickLabel',[])
% set(gca,'YTick',[-0.004,-0.002,0,0.002,0.004],'YTickLabel',[0.004,0.002,0,0.002,0.004],'fontSize',12)
axis('square');
h=legend([h1,h2],'Within-Class','Between-Class','Location','North');
rect = [0.45, 0.25, .2, .1];
set(h, 'Position', rect);

load(strcat(pre1,'ScreeningPancvsNormal.mat'))
testMGC2=testMGC;
load(strcat(pre1,'ScreeningPancvsAll.mat'))
ax=subplot('Position',F.pos5);
set(ax,'FontSize',fs-1);
% [~,ind]=sort(testMGC,'descend');
% plot(testMGC(ind),'.-','LineWidth',lw);
hold on
p1=testMGC2(:,6);
p2=testMGC(:,6);
plot(p1,p2,'.','Color',gr,'MarkerSize',5);
x=0:0.0001:1;
y=0.001;
plot(x,y*ones(length(x),1),'--','Color',glob);
plot(y*ones(length(x),1),x,'--','Color',glob);
ind=181;
plot(p1(ind),p2(ind),'.','Color',gr,'MarkerSize',15);
text(p1(ind),p2(ind),'neurogranin','VerticalAlignment','bottom','HorizontalAlignment','left','Color',gr);
hold off
set(gca,'XScale','log','YScale','log');
% xlim([1,318]);
% ylim([-0.1,0.5]);
%set(gca,'YTick',[0,0.2,0.4],'XTick',[1,100,200,300]);
% ylabel('Magnitude','FontSize',fs, ...
%     'Units', 'normalized','Position', [-0.21 0], 'HorizontalAlignment', 'left')
% xlabel('Features','FontSize',fs,...
%     'Units', 'normalized','Position', [-0.01, -0.16], 'HorizontalAlignment', 'left')
ylabel('Panc vs All','FontSize',fs, ...
    'Units', 'normalized','Position', [-0.22 0], 'HorizontalAlignment', 'left')
xlabel('Panc vs Norm','FontSize',fs,...
    'Units', 'normalized','Position', [-0.01, -0.17], 'HorizontalAlignment', 'left')
title(strcat(tits{5},'.',{' '},'Cancer Biomarker Discovery'),'FontSize',fs, ...
    'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left');
axis('square');
% semilogy(pMGC(ind));

load(strcat(rootDir,'Data/Preprocessed/BrainHippoShape.mat'));
C=LMRS;
D=squareform(pdist((Label+unifrnd(0,0.01,n,1))));
[~,per]=sort(Label);
C=C(per,per)/max(max(C));
D=D(per,per)/max(max(D));
n=size(C,1);
[A,B,RC,RD]=MGCDistTransform(C,D,'mgc');
% 
% H=eye(n)-(1/n)*ones(n,n);
% A=H*C-C/n;
% B=D*H-D/n;
% for j=1:n
%     A(j,j)=0;
%     B(j,j)=0;
% end
mcorrH=A.*B;

MH=max(max(mcorrH(2:end,2:end)))/2;
mH=min(min(mcorrH(2:end,2:end)));

cmap=zeros(2,3);
map2 = brewermap(128,'PiYG'); % brewmap
gr=map2(120,:);
pu=map2(8,:);
% cmap(1,:) = pu;
% cmap(2,:) = gr;
% map1=cmap;
% set(groot,'defaultAxesColorOrder',map1);
filename=strcat(pre1,'CorrPermDistTestTypeBrainLMRxY.mat');
load(filename);
ax=subplot('Position',F.pos3);
colormap(ax,map2);

[m,n]=size(pMLocal);
indP=optimalInd;
[J,I]=ind2sub([m,n],indP);
tmp=zeros(m,n);
tmp(J,I)=1;
tmp(testMLocal<testMGC)=0;
[k,l]=ind2sub([m,n],find(tmp==1,1,'last'));

RC=(RC<=k);
RD=(RD<=l);

ind1=reshape(RC&RD,n^2,1);
mcorrH=A.*B;
mcorrH(ind1~=1)=0;


imagesc(mcorrH');
caxis([-MH,MH]);
set(gca,'YDir','normal')
set(gca,'XTick',[1,50,100],'YTick',[1,50,100],'FontSize',fs);
title(strcat(tits{3},'.',{' '},'Joint Distance Matrix'),'FontSize',fs, ...
    'Units', 'normalized','Position', [0 1.05], 'HorizontalAlignment', 'left');
axis('square');
h=colorbar('location','eastoutside','FontSize',fs-2);
%title(h,'Magnitude')

h=suptitle(strcat({'      '},'A. Brain & Mind',{'            '},'B. Imaging Biomarkers',{'              '},'C. Screening',{'  '}));% for 1-Dimensional Simulations'));
set(h,'FontSize',15,'FontWeight','normal','Fontangle','italic','Fontname','Timesnewroman')
        
F.fname=pre2; %strcat(pre2, num2str(i));
F.wh=[8 4];
print_fig(gcf,F)

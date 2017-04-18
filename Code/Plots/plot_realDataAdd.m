function []=plot_realDataAdd(opt)
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
if nargin<1
    opt=1;
end
pre2=strcat(rootDir,'Figures/FigReal',num2str(opt));% The folder to save figures

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
    'CorrPermDistTestTypeBrainLMRxY.mat'; ...
    'CorrPermDistTestTypeMigrainxCCI.mat'};
fnames2={'BrainCP.mat'; ...
    'BrainHippoShape.mat'; ...
    'Semipar.mat'};
xlabs={ '# Activity Neighbors'; ...
        '# Shape Neighbors';...
        '# Graph Neighbors'};
ylabs={ '# Personality Neighbors'; ...
        '# Disease Neighbors';...
        '# Creativity Neighbors'};
% tits= {'A. Brain Activity vs. Personality'; ...
%     'B. Brain Shape vs. Disorder';...
%     'C. Brain Graph vs. Creativity'};
tits={'P-Value Map';'Correlation Map';'Distance Scatter Plot'};


%% loop maps
figure(1), clf, hold all
filename=strcat(pre1,fnames{2});
load(filename);
cmap2=flipud(map2);
fs=9;
cticks=[0.001, 0.01, 0.1, 0.5];
lw=1.5;

% for i=1:3
    i=opt;
    filename=strcat(pre1,fnames{i});
    load(filename);
    
    load(strcat(rootDir,'Data/Preprocessed/',fnames2{i}))
    [U,V,~]=svd(LMRS);
    X=U*V(:,1).^0.5;
    X1=X(Label==1);
    X2=X(Label==2);
    X3=X(Label==3);
    
    [~,p1,ks1] = kstest2(X1,X2);
    [~,p2,ks2] = kstest2(X1,X3);
    [~,p3,ks3] = kstest2(X2,X3);
    
    [f1,x1]=ksdensity(X1);
    [f2,x2]=ksdensity(X2);
    [f3,x3]=ksdensity(X3);
    
    hold on
    plot(x1,f1,'r')
    plot(x2,f2,'k')
    plot(x3,f3,'b')
    
    load('Semipar.mat')
    C=distMigrain(ind,ind);
    cci=cci(ind);
    nn=size(ind,1);
    [U,V]=svd(C);
    dim=10;
    X=U*V(:,1:dim).^0.5;
    X=[X ones(nn,1)];
    tn=ceil(nn/2);
    train=1:tn;
    tsn=tn+1:nn;
    rep=1;
    res=0;
    for i=1:rep
        per=randperm(nn);
        [b,~,~,~,stat] = regress(cci(per(train)),X(per(train),:));
        res=res+norm((cci(per(tsn))-X(per(tsn),:)*b)/(nn-tn))/rep;
    end
    

    [m,n]=size(pMLocal);
    
    subplot(1,3,1)
    hold on
    imagesc(log(pMLocal'));
    set(gca,'YDir','normal')
    colormap(cmap2)
    set(gca,'FontSize',fs);
    %     if i==3
    xlabel(xlabs{i},'FontSize',fs, ...
       'Units', 'normalized','Position', [-0.01, -0.14], 'HorizontalAlignment', 'left')
    ylabel(ylabs{i},'FontSize',fs, ...
        'Units', 'normalized','Position', [-0.18 0], 'HorizontalAlignment', 'left')
    title(tits{1},'FontSize',fs, ...
        'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')
    
    
    %[~,indP]=MGCScaleVerify(p2All',rep);
    indP=optimalInd;
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
    tmp(testMLocal<testMGC)=0;
    [k,l]=ind2sub([m,n],find(tmp==1,1,'last'));
    plot(k,l,'go','markerSize',6,'linewidth',3);
    plot(m,n,'.','markerSize',18,'MarkerFaceColor',glob,'Color',glob)
    xticks=[5,round(m/2)-1,m-1];
%     if i==1,  xticks(1)=3; end
    %  set(gca,'XTick',xticks,'XTickLabel',[2,round(m/2),m]); % Remove x axis ticks
%    set(gca,'YTick',[3,round(n/2)-1,n-1],'YTickLabel',[2,round(n/2),n]); % Remove x axis ticks
    xlim([2,m]);
    ylim([2,n]);
    axis('square');
    hold off
%     if i==1
%         h=colorbar('Ticks',log(cticks),'TickLabels',cticks,'location','westoutside','FontSize',fs);
%         title(h,'p-value')
%     end
% end

ax=subplot(1,3,2);
% ax=subplot('Position',[left(5), bottom(3), width, height]);
hold on
cmap2=flipud(cmap2);
set(groot,'defaultAxesColorOrder',cmap2);
% colormap(cmap2)
kmin=2;
ph=testMLocal';
%indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% ph(indPower)=2;
imagesc(ph);
plot(k,l,'go','markerSize',6,'linewidth',3)
plot(size(ph,2),size(ph,1),'.','markerSize',18,'MarkerFaceColor',glob,'Color',glob)
hold off
% set(gca,'FontSize',fs)
% set(gca,'YDir','normal')
% cmap=map1;
colormap(ax,cmap2)
% hm=ceil(max(max(ph))*100)/100;
% hm=ceil(prctile(ph(ph<1),99)*100)/100;
% caxis([0 hm])
% h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
% set(h,'FontSize',fontSize);
xlim([1 size(ph,2)]);
ylim([1 size(ph,1)]);
% set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',16);
% set(gca,'XTick',[],'YTick',[])
% pos = get(ax,'position');
axis('square');
%     xlabel(xlabs{i},'FontSize',fs, ...
%        'Units', 'normalized','Position', [-0.01, -0.14], 'HorizontalAlignment', 'left')
%     ylabel(ylabs{i},'FontSize',fs, ...
%         'Units', 'normalized','Position', [-0.16 0], 'HorizontalAlignment', 'left')
    title(tits{2},'FontSize',fs, ...
        'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')
    
ax=subplot(1,3,3);
load(strcat(rootDir,'Data/Preprocessed/',fnames2{i}))
switch i
    case 1
        C=distC;
        D=distP;     
    case 2
        C=LMRS;
        D=squareform(pdist((Label+unifrnd(0,0.01,n,1))));
    case 3
        C=distMigrain(ind,ind);
        D=squareform(pdist(cci));
        D=D(ind,ind);
end
C=C/max(max(C));
D=D/max(max(D));
% [~,per]=sort(C(:,1),'ascend');
% C=C(per,per);
% % D=squareform(pdist(y));
% D=D(per,per);
% H=eye(n)-ones(n,n)/n;
% % dcov=(H*C*H).*(H*D*H);
% for i=1:n
%     dcov(i,i)=0;
% end
% imagesc(dcov);
% set(gca,'YDir','normal')
plot(reshape(C,size(C,1)^2,1),reshape(D,size(D,1)^2,1),'k.','MarkerSize',3);
xlim([0,1]);
        ylim([0,1]);
        set(gca,'XTick',[0,1],'YTick',[0,1]);
        %     xlabel(xlabs{i},'FontSize',fs, ...
%        'Units', 'normalized','Position', [-0.01, -0.14], 'HorizontalAlignment', 'left')
%     ylabel(ylabs{i},'FontSize',fs, ...
%         'Units', 'normalized','Position', [-0.16 0], 'HorizontalAlignment', 'left')
    title(tits{3},'FontSize',fs, ...
        'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')
%         
% 
% n=42;
% % ax=subplot('Position',[left(5), bottom(3), width, height]);
% hold on
% cmap2=flipud(cmap2);
% set(groot,'defaultAxesColorOrder',cmap2);
% % colormap(cmap2)
% kmin=2;
% ph=testMLocal(kmin:n,kmin:n)';
% %indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
% % ph(indPower)=2;
% imagesc(ph);
% plot(size(ph,1)-1,size(ph,2)-1,'.','markerSize',30,'MarkerFaceColor',glob,'Color',glob)
% plot(k,l,'go','markerSize',10,'linewidth',5)
% hold off
% fontSize=13;
% set(gca,'FontSize',fontSize)
% set(gca,'YDir','normal')
% % cmap=map1;
% colormap(ax,cmap2)
% % hm=ceil(max(max(ph))*100)/100;
% % hm=ceil(prctile(ph(ph<1),99)*100)/100;
% % caxis([0 hm])
% % h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
% % set(h,'FontSize',fontSize);
% xlim([1 size(ph,1)]);
% ylim([1 size(ph,2)]);
% set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',16);
% set(gca,'XTick',[],'YTick',[])
% pos = get(ax,'position');
axis('square');
% % plot last figure
% load(strcat(pre1,'CorrBrainNoiseSummary.mat'));
% % cmap=zeros(3,3);
% % ma = [1,0,1];
% % cmap(1,:) = ma;
% % cmap(2,:) = ma;
% % map1=cmap;
% set(groot,'defaultAxesColorOrder',map1);
% 
% subplot(1,4,4)
% %scatter(x,p(:,2), 500,'k.','jitter','on', 'jitterAmount', 0.3);
% pv=p(:,1);
% [f,xi]=ksdensity(pv);
% hold on
% plot(xi,f,'.-','LineWidth',lw);
% pv=sort(pv,'ascend');
% ord=0.01*ones(length(pv),1);
% for i=2:length(pv);
%     if pv(i)-pv(i-1)<0.001
%         ord(i)=ord(i-1)+0.4;
%     end
% end
% plot(pv,ord,'.','MarkerSize',8);
% % <<<<<<< HEAD
% xlim([0,0.15]);
% ylim([-1 max(f)+1]);
% % =======
% % xlim([-0.05,0.15]);
% % ylim([-1 15]);
% % >>>>>>> 7132b2753b089edc2ca608140c119b901e31b17a
% set(gca,'FontSize',fs-2);
% set(gca,'YTick',[]); % Remove y axis ticks
% axis('square');
% hold off
% xlabel('False Positive Rate','FontSize',fs+1,...
%     'Units', 'normalized','Position', [-0.03, -0.12], 'HorizontalAlignment', 'left')
% ylabel('Density','FontSize',fs+1, ...
%     'Units', 'normalized','Position', [-0.05 0], 'HorizontalAlignment', 'left')
% % title('D. Brain Activity vs. Fake Movie','FontSize',fs+2, ...
% %     'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')
% title('D','FontSize',fs+2, ...
%     'Units', 'normalized','Position', [0 1.01], 'HorizontalAlignment', 'left')

% % F.fname=strcat(pre2, 'CORR');
% % F.wh=[3 2.5]*2;
% % print_fig(gcf,F)
% %colorbar()

h=suptitle(strcat('Brain vs Mental Properties'));% for 1-Dimensional Simulations'));
set(h,'FontSize',15,'FontWeight','normal');
        
F.fname=pre2; %strcat(pre2, num2str(i));
F.wh=[6.2 2.2];
print_fig(gcf,F)


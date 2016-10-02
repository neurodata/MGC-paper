function []=plot_schematic2
% Author: Cencheng Shen
% CorrVisualPlots()
% CorrVisualPlots(100,2)
% Used to plot figure 0 in the files

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));
pre2=strcat(rootDir,'Figures/Fig');% The folder to save figures

% if nargin<1
    n=100;
% end
% if nargin<2
    dim=1;
% end
% if nargin<3
    noise=1;
% end

total=20;
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
if dim>1
    noise=0;
end
loca=[0,1,0];
glob=[0.5,0.5,0.5];
pear=[0.2,0.2,0.2];

cmap=zeros(3,3);
gr = [0.5,0.5,0.5];
cmap(1,:) = gr;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

sz=12;
for type=1:total
    subplot(s,t,type);
    %titlechar=strcat(num2str(type),'.',{' '},CorrSimuTitle(type));
    [x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0.1);
    %[x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
    statCor=corr(x1,y1);
    C=squareform(pdist(x1));D=squareform(pdist(y1));
    statMGC=MGCSampleStat(C,D,'mcorDouble');
    A=DistCentering(C,'mcorDouble');B=DistCentering(D,'mcorDouble');
    statDcor=sum(sum(A.*B))/sqrt(sum(sum(A.*A)*sum(sum(B.*B))));
    statCor=round(statCor*1000)/1000;
    statDcor=round(statDcor*1000)/1000;
    statMGC=round(statMGC*1000)/1000;
    titlechar=strcat('\color[rgb]{.1 .1 .1}',num2str(statCor),', \color[rgb]{.5 .5 .5}',num2str(statDcor),', \color[rgb]{0 1 0}',num2str(statMGC));
    sz2=8;
    hold on
    plot(x1(:,1),y1(:,1),'k.','MarkerSize',sz2);
    if type==9;
        plot(x1(:,1),y1(:,1),'k.','MarkerSize',sz2*3);
    end
%     plot(x(:,1),y(:,1),'.','MarkerSize',sz);
    % Specify the axis limit for each type
    switch type
        case 1
            a=[-1,1];b=[-2.5,2.5];
        case 2
            a=[0,3];b=[-20,40];
        case 3
            a=[-1,1];b=[-250,250];
        case 4
            a=[-3,3];b=[-3,3];
        case 5
            a=[-1,1];b=[-2.5,2.5];
        case 6
            a=[-1,1];b=[-1,2];
        case 7
            a=[-1,1];b=[-1,2];
        case 8
            a=[-6,6];b=[-6,6];
        case 9
            a=[-0.1,1.1];b=[-2,2];
        case 10
            a=[-4,4];b=[-15,10];
        case 11
            a=[-1,1];b=[0,1.5];
        case 12
            a=[-1,1];b=[-2,2];
        case 13
            a=[-1,1];b=[-2,2];
        case 14
            a=[-2,2];b=[-2,2];
        case 15
            a=[-1,1];b=[-1.5,1.5];
        case 16
            a=[-1,1];b=[-1,1];
        case 17
            a=[-6,6];b=[-2,2];
        case 18
            a=[-2,2];b=[-2,2];
        case 19
            a=[-3,3];b=[-3,3];
        case 20
            a=[-3,3];b=[-3,3];
    end
    xlim(a);
    ylim(b);
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    ylabel(CorrSimuTitle(type),'FontSize',14);%,...
    %'Units', 'normalized','Position', [-0.01, -0.2], 'HorizontalAlignment', 'left','Interpreter','latex');
    hold off
    title(titlechar,'FontSize',16,'Interpreter','tex');
    axis('square');
%     set(gca,'box','off','ycolor','w','xcolor','w')
end
h=suptitle('Simulated Example for 20 Dependencies');
set(h,'FontSize',24,'FontWeight','normal');

F.fname=[strcat(pre2, 'SimVisual2')];
F.wh=[8 5]*2;
print_fig(gcf,F)
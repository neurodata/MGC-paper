function []=plot_simulation_visual(n,dim,noise,pre2)
% Author: Cencheng Shen
% CorrVisualPlots()
% CorrVisualPlots(100,2)
% Used to plot figure 0 in the files

%%%
fpath = mfilename('fullpath');
findex=strfind(fpath,'/');
if isempty(findex)
    findex=strfind(fpath,'\');
end
rootDir=fpath(1:findex(end-2));
p = genpath(rootDir);
gits=strfind(p,'.git');
colons=strfind(p,':');
for i=0:length(gits)-1
    endGit=find(colons>gits(end-i),1);
    p(colons(endGit-1):colons(endGit)-1)=[];
end
addpath(p);

if nargin<1
    n=100;
end
if nargin<2
    dim=1;
end
if nargin<3
    noise=1;
end
if nargin<4
    pre2='../../Figures/Fig'; % The folder to save figures
end

total=20;
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
if dim>1
    noise=0;
end

cmap=zeros(3,3);
gr = [0.5,0.5,0.5];
cmap(1,:) = gr;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

sz=12;
for type=1:total
    subplot(s,t,type);
    titlechar=CorrSimuTitle(type);
    [x, y]=CorrSampleGenerator(type,n,dim,1, noise);
    [x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
    if type==19
        sz2=24; % Enlarge the point size for discrete distribution, i.e., uncorrelation binomial
    else
        sz2=8;
    end
    hold on
    plot(x1(:,1),y1(:,1),'k.','MarkerSize',sz2);
    plot(x(:,1),y(:,1),'.','MarkerSize',sz);
    % Specify the axis limit for each type
    switch type
        case 1
            a=[-1,1];b=[-2.5,2.5];
        case 2
            a=[-1,1];b=[-250,250];
        case 3
            a=[-1,1];b=[-2.5,2.5];
        case 4
            a=[0,3];b=[-20,40];
        case 5
            a=[-3,3];b=[-3,3];
        case 6
            a=[-1,1];b=[-1,2];
        case 7
            a=[-1,1];b=[-1,2];
        case 8
            a=[-1,1];b=[-1.5,1.5];
        case 9
            a=[-1,1];b=[0,1.5];
        case 10
            a=[-4,4];b=[-15,10];
        case 11
            a=[-1,1];b=[-1,1];
        case 12
            a=[-6,6];b=[-2,2];
        case 13
            a=[-6,6];b=[-6,6];
        case 14
            a=[-2,2];b=[-2,2];
        case 15
            a=[-2,2];b=[-2,2];
        case 16
            a=[-1,1];b=[-2,2];
        case 17
            a=[-1,1];b=[-2,2];
        case 18
            a=[-3,3];b=[-3,3];
        case 19
            a=[-0.1,1.1];b=[-2,2];
        case 20
            a=[-3,3];b=[-3,3];
    end
    xlim(a);
    ylim(b);
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    hold off
    title(titlechar,'FontSize',14);
end
h=suptitle('Visualization for 20 Simulated Dependencies');
set(h,'FontSize',24,'FontWeight','normal');

F.fname=[strcat(pre2, 'SimVisual')];
F.wh=[8 4]*2;
print_fig(gcf,F)

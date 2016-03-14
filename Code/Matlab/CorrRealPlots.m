function []=CorrRealPlots(pre1,pre2)
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

if nargin<1
    pre1='../../Data/'; % The folder to locate data
end
if nargin<2
    pre2='../../Figures/FigReal'; % The folder to save figures
end
cmap=zeros(3,3);
gr = [0,1,0];
ma = [1,0,1];
cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = ma;
cmap(3,:) = cy;
map1=cmap;
map2 = brewermap(128,'GnBu'); % brewmap

% Plot heatmap
total=4;
set(groot,'defaultAxesColorOrder',map1);
[filename, titleStr]=CorrRealDataName(1);
filename=strcat(pre1,filename);
load(filename);
n=size(p1All,1);
kmin=2;
xa=kmin:n;
pp1=p1All(xa,4);
[filename, titleStr]=CorrRealDataName(2);
filename=strcat(pre1,filename);
load(filename);
n=size(p1All,1);
pp2=p1All(xa,4);
figure
plot(xa,pp1,'.-',xa,pp2,'.-','LineWidth',2)
xlabel('Neighborhood Choice of X','FontSize',16);
xlim([1 n]);
ylim([0 0.2]);
ylabel('P-Value','FontSize',16);

title('Local Tests P-value for Brain Shape vs Disorder','FontSize',13);
h=legend('Left Brain Shape vs Disorder','Right Brain Shape vs Disorder');
set(h,'FontSize',13);

F.fname=strcat(pre2, num2str(1));
F.wh=[3 2.5]*2;
print_fig(gcf,F)

for i=3:total
    [filename, titleStr]=CorrRealDataName(i);
    filename=strcat(pre1,filename);
    load(filename);
    n=size(p1All,1);
    figure
    kmin=2;
    imagesc(p1All(kmin:n,kmin:n)');
    set(gca,'YDir','normal')
    colormap(flipud(map2))
    caxis([0.01 0.1])
    colorbar
    xlabel('Neighborhood Choice of X','FontSize',16);
    ylabel('Neighborhood Choice of Y','FontSize',16);
    
    % Figure title/labels
    titleStr = strcat('Local Tests P-value for ', titleStr);
    title(titleStr,'FontSize',13);
    
    F.fname=strcat(pre2, num2str(i));
    F.wh=[3 2.5]*2;
    print_fig(gcf,F)
end

function [str, title]=CorrRealDataName(i)
str='CorrPermDistTestType';
switch i
%     case 1
%         str=strcat(str,'BrainCxP.mat');
%         title=' Connectome vs Personality';
    case 1
        str=strcat(str,'BrainLMLxY.mat');
        title=' Left Brain Shape vs Disorder';
    case 2
        str=strcat(str,'BrainLMRxY.mat');
        title=' Right Brain Shape vs Disorder';
    case 3
        str=strcat(str,'MigrainxCCI.mat');
        title=' Migrain vs CCI';
    case 4
        str=strcat(str,'M2gxCCI.mat');
        title=' M2g vs CCI';
end
function []=CorrRealPlots(pre1,pre2)
% Used to plot the heatmap of real data used in tex. Run like
% CorrRealPlots()

if nargin<1
    pre1='../../Data/'; % The folder to locate data
end
if nargin<2
    pre2='../../Figures/FigReal'; % The folder to save figures
end

% Plot heatmap
total=3;
for i=1:total
    [filename, titleStr]=CorrRealDataName(i);
    filename=strcat(pre1,filename);
    load(filename);
    figure
    K=n;
    kmin=1;
    if n>50
        c=2;
        K=ceil(K/2);
    else
        c=1;
        kmin=2;
    end
    xaxis=kmin:K;
    yaxis=kmin:K;
    [X,Y]=meshgrid(c*xaxis,c*yaxis);
    ph=p1All(c*xaxis,c*yaxis)';
    surf(X,Y,ph);
    view(2)
    colormap(flipud(colormap))
    caxis([0.01 0.1])
    colorbar
    xlabel('Neighborhood Choice of X','FontSize',16);
    ylabel('Neighborhood Choice of Y','FontSize',16);
    xlim([c*kmin,c*K]);
    ylim([c*kmin,c*K]);
    
    % Figure title/labels
    titleStr = strcat('P-values of All Local Tests for ', titleStr);
    title(titleStr,'FontSize',13);
    
    F.fname=strcat(pre2, num2str(i));
    F.wh=[3 2.5]*2;
    print_fig(gcf,F)
end

function [str, title]=CorrRealDataName(i)
str='CorrPermDistTestType';
switch i
    case 1
        str=strcat(str,'BrainCxP.mat');
        title=' Connectome vs Personality';
    case 2
        str=strcat(str,'BrainLMLxY.mat');
        title=' Left Brain Shape vs Disorder';
    case 3
        str=strcat(str,'BrainLMRxY.mat');
        title=' Right Brain Shape vs Disorder';
end
function []=CorrVisualPlots(n,dim,noise,pre2)
% Author: Cencheng Shen
% CorrVisualPlots()
% Used to plot figure 0 in the files
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
    pre2='../../Figures/JovoFig'; % The folder to save figures
    %pre2='News_1/Fig';
end

total=20;
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for type=1:total
    subplot(s,t,type);
    titlechar=CorrSimuTitle(type);
    [x, y]=CorrSampleGenerator(type,n,dim,1, noise);
    [x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency 
    if type==19
        sz=20; % Enlarge the point size for discrete distribution, i.e., uncorrelation binomial
    else
        sz=2;
    end
        hold on
        plot(x1,y1,'r.','MarkerSize',sz);
        plot(x,y,'.');
        hold off
    title(titlechar);
end
h=suptitle('Visualization for 20 Simulated Dependencies');
set(h,'FontSize',20,'FontWeight','normal');

F.fname=[strcat(pre2, '0')];
F.wh=[8 4]*2;
print_fig(gcf,F)

function []=CorrVisualPlots(n,dim,noise)
% Author: Cencheng Shen
% n=100;dim=1;noise=1;
% CorrVisualPlots(n,dim,noise)
% Used to plot figure 0 in the files
total=20;

figure
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
suptitle('Visualization for 20 Simulated Dependencies')

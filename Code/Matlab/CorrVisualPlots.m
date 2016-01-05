function []=CorrVisualPlots(n,dim,noise)
% Author: Cencheng Shen
% n=1000;dim=1;noise=0;
% CorrVisualPlots(n,dim,noise)
% Used to plot figure 0 in the files
total=20;
%noiseC=noise*ones(total,1);
%noiseC(12)=0;noiseC(13)=0;noiseC(15)=0;noiseC(17)=0;noiseC(18)=0;

figure
s=4;
t=5;
for type=1:total
    subplot(s,t,type);
    titlechar=CorrSimuTitle(type);
    [x, y]=CorrSampleGenerator(type,n,dim,1, noise);
    [x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0);
    if type==18
        sz=20;
    else
        sz=2;
    end
        hold on
        plot(x1,y1,'r.','MarkerSize',sz);
        plot(x,y,'.');
        hold off
%     end
    title(titlechar);
end
suptitle('Visualization for 20 Simulated Dependencies')

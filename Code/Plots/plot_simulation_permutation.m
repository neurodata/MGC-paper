function []=plot_simulation_permutation

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data
pre2=strcat(rootDir,'Figures/Fig');% The folder to save figures

%% Set colors
%map2 = brewermap(128,'PRGn'); % brewmap
loca=[0,1,0];
glob= [1,0,1];
HHG   = [0.5,0.5,0.5];
mk=30;
% map1=zeros(7,3);
% gr = [0,1,0];
% ma = [1,0,1];
% cy = [0,1,1];
% cmap(1,:) = gr;
% cmap(2,:) = ma;
% cmap(3,:) = cy;
% cmap(7,:) = [0 0 0];
% dcorr = cmap(1,:);
% mcorr = cmap(2,:);
% mante = cmap(3,:);
% HHG   = [0.5,0.5,0.5];

ls{1}='-';
ls{2}='--';
ls{3}=':';
ls{4}='--';
ls{1}='.';
ls{2}='.';
ls{3}='.';
ls{4}='.';

%
% map1(1,:)=dcorr;
% map1(2,:)=mcorr; 
% map1(3,:)=mcorr; 
% map1(4,:)=HHG; 
% set(groot,'defaultAxesColorOrder',map1);

figNumber='1DPerm';
filename=strcat(pre1,'CorrSimPermScale1-20Dim1');
load(filename)
figure
x=1:20;
a=2;
if a==1
    ind=[1,2,3,7];
else
    ind=[4,5,6,7];
end
p1=powerP(ind,:);
p1=p1-repmat(p1(1,:),4,1);
hold on
plot(x,p1(3,:),ls{1},'Color',loca,'LineWidth',4,'MarkerSize',mk);
plot(x,p1(2,:),ls{1},'Color',glob,'LineWidth',4,'MarkerSize',mk);
plot(x,p1(4,:),ls{4},'Color',HHG,'LineWidth',4,'MarkerSize',mk);
hold off
legend('True MGC','Mcorr','HHG','Location','NorthWest');
set(gca,'FontSize',24);
legend boxoff
xlabel('Function Type','FontSize',24);
ylabel('Power Difference','FontSize',24);
xlim([1,20]);
ylim([-0.6,1]);
set(gca,'XTick',[1,10,20]); % Remove x axis ticks
    set(gca,'YTick',[-0.5,0,0.5,1]); % Remove x axis ticks
title('1-Dimensional Simulations','FontSize',24);
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

figNumber='HDPerm';
filename=strcat(pre1,'CorrSimPermScale1-20Dim2');
load(filename)
figure
x=1:20;
a=2;
if a==1
    ind=[1,2,3,7];
else
    ind=[4,5,6,7];
end
p1=powerP(ind,:);
p1=p1-repmat(p1(1,:),4,1);
hold on
plot(x,p1(3,:),ls{1},'Color',loca,'LineWidth',4,'MarkerSize',mk);
plot(x,p1(2,:),ls{1},'Color',glob,'LineWidth',4,'MarkerSize',mk);
plot(x,p1(4,:),ls{4},'Color',HHG,'LineWidth',4,'MarkerSize',mk);
hold off
%legend('True MGC','Mcorr','HHG','Location','NorthWest');
set(gca,'FontSize',24);
xlabel('Function Type','FontSize',24);
%ylabel('Power Difference','FontSize',24);
xlim([1,20]);
ylim([-0.6,1]);
set(gca,'XTick',[1,10,20]); % Remove x axis ticks
set(gca,'YTick',[-0.5,0,0.5,1]); % Remove x axis ticks
title('High-Dimensional Simulations','FontSize',24);
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)
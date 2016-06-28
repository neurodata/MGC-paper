function []=plot_simulation_permutation2

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

ls{1}='-';
ls{2}='-';
ls{3}='--';

% set(groot,'defaultAxesColorOrder',map1);

%% 1-d
figNumber='1DPerm';
filename=strcat(pre1,'CorrSimPermScale1-20Dim1');
load(filename)

x=1:20;
a=2;
if a==1
    ind=[1,2,3,7];
else
    ind=[4,5,6,7];
end
p1=powerP(ind,:);
p1=p1-repmat(p1(1,:),4,1);
pp1=p1(2,:);
pp2=p1(3,:);
pp3=p1(4,:);
[f1,xi1]=ksdensity(pp1,'support',[-1,1]);
[f2,xi2]=ksdensity(pp2,'support',[-1,1]);
[f3,xi3]=ksdensity(pp3,'support',[-1,1]);

figure(1)%, clf
hold on
h3=plot(x,pp3,ls{3},'Color',HHG,'LineWidth',4);
h2=plot(x,pp2,ls{2},'Color',glob,'LineWidth',4);
h1=plot(x,pp1,ls{1},'Color',loca,'LineWidth',4);
hold off
h=legend([h1 h2 h3],'True MGC','Mcorr','HHG','Location','NorthWest');
set(gca,'FontSize',24);

legend boxoff
% xlabel('Function Type','FontSize',24);
%ylabel('Power Difference','FontSize',24);
xlim([1,20]);
ylim([-0.6,1]);
set(gca,'XTick',[1,10,20]); % Remove x axis ticks
    set(gca,'YTick',[-0.5,0,0.5,1]); % Remove x axis ticks
title('1-Dimensional Simulations','FontSize',24);
grid on
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
% print_fig(gcf,F)

%% high dim
figNumber='HDPerm';
filename=strcat(pre1,'CorrSimPermScale1-20Dim2');
load(filename)
figure(2), clf
x=1:20;
a=2;
if a==1
    ind=[1,2,3,7];
else
    ind=[4,5,6,7];
end
p1=powerP(ind,:);
p1=p1-repmat(p1(1,:),4,1);
pp1=p1(2,:);
pp2=p2(2,:);
pp3=p3(2,:);
[f1,xi1]=ksdensity(p1(2,:),'support',[-1,1]);
[f2,xi2]=ksdensity(p1(3,:),'support',[-1,1]);
[f3,xi3]=ksdensity(p1(4,:),'support',[-1,1]);

[~, so]=sort(p1(3,:),'descend');

hold on
plot(x,p1(3,so),ls{4},'Color',loca,'LineWidth',4,'MarkerSize',mk);
plot(x,p1(2,so),ls{5},'Color',glob,'LineWidth',4,'MarkerSize',mk);
plot(x,p1(4,so),ls{6},'Color',HHG,'LineWidth',4,'MarkerSize',mk);
plot([0 20], [0 0],'color',0.5*[1 1 1])
hold off
%legend('True MGC','Mcorr','HHG','Location','NorthWest');
set(gca,'FontSize',24);
% xlabel('Function Type','FontSize',24);
%ylabel('Power Difference','FontSize',24);
xlim([1,20]);
ylim([-0.6,1]);
set(gca,'XTick',[1,10,20]); % Remove x axis ticks
set(gca,'YTick',[-0.5,0,0.5,1]); % Remove x axis ticks
title('High-Dimensional Simulations','FontSize',24);
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
% print_fig(gcf,F)

%%

% 1 estimated mgc
% 2 mcorr
% 3 true mgc
% 4 hhg

figure(4), clf
hold on
p2=powerP(ind,:);
bw=0.3;
jit=inf;
sp=-0.3;

% True MGC
[foo,bar]=ksdensity(p2(3,:),'support',[0 1],'bandwidth',bw);
h(1)=plot(bar,foo,ls{6},'Color','k','LineWidth',4,'MarkerSize',mk);
plot(p2(3,:),4*sp+zeros(20,1)+rand(20,1)/jit,ls{5},'Color','k','LineWidth',4,'MarkerSize',mk);

% estimated MGC
[foo,bar]=ksdensity(p2(1,:),'support',[0 1],'bandwidth',bw);
h(2)=plot(bar,foo,ls{6},'Color',loca,'LineWidth',4,'MarkerSize',mk);
plot(p2(1,:),3*sp+zeros(20,1)+rand(20,1)/jit,ls{5},'Color',loca,'LineWidth',4,'MarkerSize',mk);

% Mcorr
[foo,bar]=ksdensity(p2(2,:),'support',[0 1],'bandwidth',bw);
h(3)=plot(bar,foo,ls{6},'Color',glob,'LineWidth',4,'MarkerSize',mk);
plot(p2(2,:),2*sp+zeros(20,1)+rand(20,1)/jit,ls{5},'Color',glob,'LineWidth',4,'MarkerSize',mk);

% HHG
[foo,bar]=ksdensity(p2(4,:),'support',[0 1],'bandwidth',bw);
h(4)=plot(bar,foo,ls{6},'Color',HHG,'LineWidth',4,'MarkerSize',mk);
plot(p2(4,:),1*sp+zeros(20,1)+rand(20,1)/jit,ls{5},'Color',HHG,'LineWidth',4,'MarkerSize',mk);

plot([0 1], [0, 0],'--k')
hold off
% legend(h,'True MGC','Estimated MGC','Mcorr','HHG','Location','NorthEast');
set(gca,'FontSize',24);
legend boxoff
ylabel('Probability','FontSize',24);
xlabel('Power','FontSize',24);
% xlim([1,20]);
% ylim([-0.6,1]);
set(gca,'XTick',[0,0.5,1]); % Remove x axis ticks
set(gca,'YTick',[0]); % Remove x axis ticks
title('High-Dimensional Simulations','FontSize',24);
% grid on
axis('tight')
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%

% 1 estimated mgc
% 2 mcorr
% 3 true mgc
% 4 hhg

figure(5), clf
hold on
p2=powerP(ind,:);
p2=p2-repmat(p2(3,:),4,1);
bw=0.3;

[foo,bar]=ksdensity(p2(1,:),'support',[-1 1],'bandwidth',bw);
plot(bar,foo,ls{6},'Color',loca,'LineWidth',4,'MarkerSize',mk);

[foo,bar]=ksdensity(p2(2,:),'support',[-1 1],'bandwidth',bw);
plot(bar,foo,ls{6},'Color',glob,'LineWidth',4,'MarkerSize',mk);

% [foo,bar]=ksdensity(p2(3,:),'support',[-1 1],'bandwidth',bw);
% plot(bar,foo,ls{6},'Color','k','LineWidth',4,'MarkerSize',mk);

[foo,bar]=ksdensity(p2(4,:),'support',[-1 1],'bandwidth',bw);
plot(bar,foo,ls{6},'Color',HHG,'LineWidth',4,'MarkerSize',mk);
hold off
legend('Estimated MGC','Mcorr','True MGC','HHG','Location','NorthEast');
set(gca,'FontSize',14);
legend boxoff
ylabel('Probability','FontSize',24);
xlabel('Power','FontSize',24);
% xlim([1,20]);
% ylim([-0.6,1]);
set(gca,'XTick',[-1:0.5:1]); % Remove x axis ticks
% set(gca,'YTick',[-0.5,0,0.5,1]); % Remove x axis ticks
title('High-Dimensional Simulations','FontSize',24);
% grid on
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
% print_fig(gcf,F)


%%

figure(6), clf
hold on
p2=powerP(ind,:);
plot(-1+ones(20,1)+randn(20,1)/10,p2(3,:),ls{5},'Color','k','LineWidth',4,'MarkerSize',mk);
plot(ones(20,1)+randn(20,1)/10,p2(1,:),ls{5},'Color',loca,'LineWidth',4,'MarkerSize',mk);
plot(1+ones(20,1)+randn(20,1)/10,p2(2,:),ls{5},'Color',glob,'LineWidth',4,'MarkerSize',mk);
plot(2+ones(20,1)+randn(20,1)/10,p2(4,:),ls{5},'Color',HHG,'LineWidth',4,'MarkerSize',mk);
hold off
% legend('True MGC','Estimated MGC','Mcorr','HHG','Location','NorthWest');
set(gca,'FontSize',24);
legend boxoff
% xlabel('Function Type','FontSize',24);
ylabel('Power','FontSize',24);
% xlim([1,20]);
% ylim([-0.6,1]);
set(gca,'XTick',[0,1,2,3],'XtickLabel',[{'True MGC'}; {'Estimated MGC'}; {'Mcorr'}; {'HHG'}]); % Remove x axis ticks
set(gca,'YTick',[0,0.5,1]); % Remove x axis ticks
ax=gca;
ax.XTickLabelRotation=45;
title('High-Dimensional Simulations','FontSize',24);
grid on
xlim([-0.3, 3.3])
F.fname=strcat(pre2, figNumber, 'A');
F.wh=[3 2.5]*2;
print_fig(gcf,F)
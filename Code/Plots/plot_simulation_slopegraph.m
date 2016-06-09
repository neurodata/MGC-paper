function []=plot_simulation_slopegraph(pre1,pre2)
% Used to plot figure 1-8 used in tex. Run like

% total is usually 20.
% pre1 specifies the location to load data.
% pre2 specifies the location to save pictures.

if nargin<1
    pre1='../../Data/Results/'; % The folder to locate data
end
if nargin<2
    pre2='../../Draft/Figures/Fig'; % The folder to save figures
end
total=20;

%% Set colors
map1=zeros(7,3);
gr = [0,1,0];
ma = [1,0,1];
cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = ma;
cmap(3,:) = cy;
cmap(7,:) = [0 0 0];
dcorr = cmap(1,:);
mcorr = cmap(2,:);
mante = cmap(3,:);
HHG   = [0.5,0.5,0.5];
map1(1,:)=dcorr; map1(5,:)=dcorr; % The color for MGC{dcorr} and global dcorr.
map1(2,:)=mcorr; map1(6,:)=mcorr; % The color for MGC{mcorr} and global mcorr.
map1(3,:)=mante; map1(7,:)=mante; % The color for MGC{Mantel} and global Mantel.
map1(4,:)=HHG; % The color for HHG
set(groot,'defaultAxesColorOrder',map1);
strL=['Global';'Local '];

%% %performance profile
figNumber='1DSlope';
figure
set(groot,'defaultAxesColorOrder',map1)
a=0;b=1;
AUC=zeros(7,1);
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    mm=1;%max([mean(power1), mean(power2),mean(power3),mean(power4),mean(power5),mean(power6),mean(power7)]);
    AUC(1)=AUC(1)+mean(power1)/mm;
    AUC(2)=AUC(2)+mean(power2)/mm;
    AUC(3)=AUC(3)+mean(power3)/mm;
    AUC(4)=AUC(4)+mean(power4)/mm;
    AUC(5)=AUC(5)+mean(power5)/mm;
    AUC(6)=AUC(6)+mean(power6)/mm;
    AUC(7)=AUC(7)+mean(power7)/mm;
end
AUC=AUC./total;
x=1:2;

sumP=zeros(4,2); sumP(1,1)=AUC(4);sumP(1,2)=AUC(1);sumP(2,1)=AUC(5);sumP(2,2)=AUC(2);sumP(3,1)=AUC(6);sumP(3,2)=AUC(3);sumP(4,1)=AUC(4);sumP(4,2)=AUC(4);
plot(x,sumP(1,:),'.-',x,sumP(2,:),'.-',x,sumP(3,:),'.--',x,sumP(4,:),'.--','LineWidth',3);
legend('dcorr','mcorr','Mantel','HHG','Location','NorthWest');
legend boxoff
% ln = findobj('type','line');
% set(ln,'marker','.','markers',14,'markerfa','w')
% text(x(:,1),x(:,2),'1')
%ylabel('power');
ylim([0,1]);
yTickN=[floor(sumP(1,1)*100)/100,1];
set(gca,'XTickLabel',strL,'XTick',1:2,'FontSize',16,'YTick',yTickN);
yTickN=[floor(sumP(1,2)*100)/100,1];
axes('xlim', [1 2],'ylim', [0 1], 'color', 'none', 'YAxisLocation', 'right','XTick',[],'FontSize',16);
set(gca,'YTick',yTickN);
%axes( 'YTick',yTickN,'YAxisLocation', 'right')
title('Mean Powers for 1-Dimensional Simulations');
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%%performance profile AUC for n
figNumber='HDSlope';
figure
set(groot,'defaultAxesColorOrder',map1)
a=0;b=1;
AUC=zeros(7,1);
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    mm=1;%max([mean(power1), mean(power2),mean(power3),mean(power4),mean(power5),mean(power6),mean(power7)]);
    AUC(1)=AUC(1)+mean(power1)/mm;
    AUC(2)=AUC(2)+mean(power2)/mm;
    AUC(3)=AUC(3)+mean(power3)/mm;
    AUC(4)=AUC(4)+mean(power4)/mm;
    AUC(5)=AUC(5)+mean(power5)/mm;
    AUC(6)=AUC(6)+mean(power6)/mm;
    AUC(7)=AUC(7)+mean(power7)/mm;
end
AUC=AUC./total;
x=1:2;

sumP=zeros(4,2); sumP(1,1)=AUC(4);sumP(1,2)=AUC(1);sumP(2,1)=AUC(5);sumP(2,2)=AUC(2);sumP(3,1)=AUC(6);sumP(3,2)=AUC(3);sumP(4,1)=AUC(4);sumP(4,2)=AUC(4);
plot(x,sumP(1,:),'.-',x,sumP(2,:),'.-',x,sumP(3,:),'.--',x,sumP(4,:),'.--','LineWidth',3);
legend('dcorr','mcorr','Mantel','HHG','Location','NorthWest');
legend boxoff
% ln = findobj('type','line');
% set(ln,'marker','.','markers',14,'markerfa','w')
% text(x(:,1),x(:,2),'1')
%ylabel('power');
ylim([0,1]);
yTickN=[floor(sumP(3,1)*100)/100,floor(sumP(2,1)*100)/100,1];
set(gca,'XTickLabel',strL,'XTick',1:2,'FontSize',16,'YTick',yTickN);
yTickN=[floor(sumP(3,2)*100)/100,floor(sumP(2,2)*100)/100,1];
axes('xlim', [1 2],'ylim', [0 1], 'color', 'none', 'YAxisLocation', 'right','XTick',[],'FontSize',16);
set(gca,'YTick',yTickN);
%axes( 'YTick',yTickN,'YAxisLocation', 'right')
title('Mean Powers for High-Dimensional Simulations');
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

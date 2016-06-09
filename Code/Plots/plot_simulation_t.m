function []=plot_simulation_t(pre1,pre2)
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
select=1;
if select==1
    map1(1,:)=dcorr; map1(4,:)=dcorr; % The color for MGC{dcorr} and global dcorr.
    map1(2,:)=mante; map1(5,:)=mante; % The color for MGC{Mantel} and global Mantel.
    map1(3,:)=HHG; % The color for HHG
end

set(groot,'defaultAxesColorOrder',map1);
lim=20;
%% %performance profile
set(groot,'defaultAxesColorOrder',map1)
a=0;b=1;
AUC1=zeros(40,lim);
AUC2=zeros(40,lim);
AUC3=zeros(40,lim);
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    mm=1;%max([mean(power1), mean(power2),mean(power3),mean(power4),mean(power5),mean(power6),mean(power7)]);
    AUC1(j,:)=power1-power7;
    AUC1(20+j,:)=power4-power7;
    AUC2(j,:)=power2-power7;
    AUC2(20+j,:)=power5-power7;
    AUC3(j,:)=power3-power7;
    AUC3(20+j,:)=power6-power7;
end

figure
hold on
x=1:lim;
for j=1:20
    plot(x,AUC1(j,:),'.:','LineWidth',1,'Color',dcorr/2);
end
for j=21:40
    plot(x,AUC1(j,:),'.-','LineWidth',1,'Color',dcorr);
end
hold off
figNumber='1DDcorr';
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

figure
hold on
x=1:lim;
for j=1:20
    plot(x,AUC2(j,:),'.:','LineWidth',1,'Color',mcorr/2);
end
for j=21:40
    plot(x,AUC2(j,:),'.-','LineWidth',1,'Color',mcorr);
end
hold off
figNumber='1DMcorr';
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

figure
hold on
x=1:lim;
for j=1:20
    plot(x,AUC3(j,:),'.:','LineWidth',1,'Color',mante/2);
end
for j=21:40
    plot(x,AUC3(j,:),'.-','LineWidth',1,'Color',mante);
end
hold off
figNumber='1DMantel';
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)



AUC1=zeros(40,lim);
AUC2=zeros(40,lim);
AUC3=zeros(40,lim);
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    if length(dimRange)<20
        continue;
    end
    mm=1;%max([mean(power1), mean(power2),mean(power3),mean(power4),mean(power5),mean(power6),mean(power7)]);
    AUC1(j,:)=power1-power7;
    AUC1(20+j,:)=power4-power7;
    AUC2(j,:)=power2-power7;
    AUC2(20+j,:)=power5-power7;
    AUC3(j,:)=power3-power7;
    AUC3(20+j,:)=power6-power7;
end

figure
hold on
x=1:lim;
for j=1:20
    plot(x,AUC1(j,:),'.:','LineWidth',1,'Color',dcorr/2);
end
for j=21:40
    plot(x,AUC1(j,:),'.-','LineWidth',1,'Color',dcorr);
end
hold off
figNumber='HDDcorr';
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

figure
hold on
x=1:lim;
for j=1:20
    plot(x,AUC2(j,:),'.:','LineWidth',1,'Color',mcorr/2);
end
for j=21:40
    plot(x,AUC2(j,:),'.-','LineWidth',1,'Color',mcorr);
end
hold off
figNumber='HDMcorr';
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

figure
hold on
x=1:lim;
for j=1:20
    plot(x,AUC3(j,:),'.:','LineWidth',1,'Color',mante/2);
end
for j=21:40
    plot(x,AUC3(j,:),'.-','LineWidth',1,'Color',mante);
end
hold off
figNumber='HDMantel';
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)
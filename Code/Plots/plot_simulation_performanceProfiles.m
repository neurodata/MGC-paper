function []=plot_simulation_powers(pre1,pre2)
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

%% %performance profile
figNumber='1DPP';
figure
set(groot,'defaultAxesColorOrder',map1)
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
profile=zeros(7,length(xaxis));
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    thres=0.8;
    ind=[find(power2>=thres,1) find(power4>=thres,1) find(power5>=thres,1) find(power6>=thres,1) find(power7>=thres,1) lim];
    pos=min(ind);
    power=[power1(pos), power2(pos), power3(pos), power4(pos),power5(pos),power6(pos),power7(pos)];
    pmax=max(power);
    for k=1:length(power)
        tmp=(pmax-power(k));
        tmpInd=ceil(tmp/interval)+1;
        profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
    end
end
profile=profile./total;
sumP=ceil(mean(profile,2)*1000)/1000;
plot(xaxis,profile(1,:),'.-',xaxis,profile(2,:),'.-',xaxis,profile(3,:),'.-',xaxis,profile(7,:),'.--',xaxis,profile(4,:),'.:',xaxis, profile(5,:),'.:',xaxis,profile(6,:),'.:','LineWidth',2);
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
set(gca,'FontSize',13);
legend(strcat('MGC\{dcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{mcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{Mantel\}, AUC=', num2str(sumP(3))),strcat('HHG, AUC=', num2str(sumP(7))),strcat('dcorr, AUC=', num2str(sumP(4))),strcat('mcorr, AUC=', num2str(sumP(5))),strcat('Mantel, AUC=', num2str(sumP(6))),'Location','SouthEast');
legend boxoff
xlabel('Difference with the Best Method','FontSize',15);
ylabel('Relative Performance','FontSize',15);
ylim([0 1]);
titleStr = strcat('Performance Profiles for 1-Dimensional Settings');
title(titleStr,'FontSize',15);
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%%performance profile AUC for n
figNumber='1DPPAUC';
figure
set(groot,'defaultAxesColorOrder',map1)
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
limN=10;
sumP=zeros(7,limN);
%load data
for ll=1:limN
    profile=zeros(7,length(xaxis));
    thres=ll/limN;
    for j=1:total
        filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
        load(filename)
        ind=[find(power2>=thres,1) find(power4>=thres,1) find(power5>=thres,1) find(power6>=thres,1) find(power7>=thres,1) lim];
        pos=min(ind);
        power=[power1(pos), power2(pos), power3(pos), power4(pos),power5(pos),power6(pos),power7(pos)];
        pmax=max(power);
        for k=1:length(power)
            tmp=(pmax-power(k));
            tmpInd=ceil(tmp/interval)+1;
            profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
        end
    end
    profile=profile./total;
    sumP(:,ll)=ceil(mean(profile,2)*1000)/1000;
end
xaxis=1/limN:1/limN:1;
plot(xaxis,sumP(1,:),'.-',xaxis,sumP(2,:),'.-',xaxis,sumP(3,:),'.-',xaxis,sumP(7,:),'.--',xaxis, sumP(4,:),'.:',xaxis,sumP(5,:),'.:',xaxis, sumP(6,:),'.:','LineWidth',2);
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
legend('MGC\{dcorr\}','MGC\{mcorr\}','MGC\{Mantel\}','HHG','dcorr','mcorr','Mantel','Location','SouthWest');
set(gca,'FontSize',13);
legend boxoff
xlabel('Threshold of Power','FontSize',15);
ylabel('Area Under Curve','FontSize',15);
xlim([1/limN 1]);
ylim([0 1]);
titleStr = strcat('AUC of Performance Profiles for 1-Dimensional Settings');
title(titleStr,'FontSize',14);
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%%performance profile
figNumber='HDPP';
figure
set(groot,'defaultAxesColorOrder',map1)
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
profile=zeros(7,length(xaxis));
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    thres=0.5;
    ind=[find(power2>=thres,1,'last') find(power4>=thres,1,'last') find(power5>=thres,1,'last') find(power6>=thres,1,'last') find(power7>=thres,1,'last') 1];
    pos=max(ind);
    power=[power1(pos), power2(pos), power3(pos), power4(pos),power5(pos),power6(pos),power7(pos)];
    pmax=max(power);
    for k=1:length(power)
        tmp=(pmax-power(k));
        tmpInd=ceil(tmp/interval)+1;
        profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
    end
end
profile=profile./total;
sumP=ceil(mean(profile,2)*1000)/1000;
plot(xaxis,profile(1,:),'.-',xaxis,profile(2,:),'.-',xaxis,profile(3,:),'.-',xaxis,profile(7,:),'.--',xaxis,profile(4,:),'.:',xaxis, profile(5,:),'.:',xaxis,profile(6,:),'.:','LineWidth',2);
legend(strcat('MGC\{dcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{mcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{Mantel\}, AUC=', num2str(sumP(3))),strcat('HHG, AUC=', num2str(sumP(7))),strcat('dcorr, AUC=', num2str(sumP(4))),strcat('mcorr, AUC=', num2str(sumP(5))),strcat('Mantel, AUC=', num2str(sumP(6))),'Location','SouthEast');
set(gca,'FontSize',13);
legend boxoff
xlabel('Difference with the Best Method','FontSize',15);
ylabel('Relative Performance','FontSize',15);
ylim([0 1]);
titleStr = strcat('Performance Profiles for High-Dimensional Settings');
title(titleStr,'FontSize',15);
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%%performance profile
figNumber='HDPPAUC';
figure
set(groot,'defaultAxesColorOrder',map1)
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
limN=10;
sumP=zeros(7,limN);
%load data
for ll=1:limN
    profile=zeros(7,length(xaxis));
    thres=ll/limN;
    for j=1:total
        filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
        load(filename)
        ind=[find(power2>=thres,1,'last') find(power4>=thres,1,'last') find(power5>=thres,1,'last') find(power6>=thres,1,'last') find(power7>=thres,1,'last') 1];
        pos=max(ind);
        power=[power1(pos), power2(pos), power3(pos), power4(pos),power5(pos),power6(pos),power7(pos)];
        pmax=max(power);
        for k=1:length(power)
            tmp=(pmax-power(k));
            tmpInd=ceil(tmp/interval)+1;
            profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
        end
    end
    profile=profile./total;
    sumP(:,ll)=ceil(mean(profile,2)*1000)/1000;
end
xaxis=1/limN:1/limN:1;
plot(xaxis,sumP(1,:),'.-',xaxis,sumP(2,:),'.-',xaxis,sumP(3,:),'.-',xaxis,sumP(7,:),'.--',xaxis, sumP(4,:),'.:',xaxis,sumP(5,:),'.:',xaxis, sumP(6,:),'.:','LineWidth',2);
legend('MGC\{dcorr\}','MGC\{mcorr\}','MGC\{Mantel\}','HHG','dcorr','mcorr','Mantel','Location','SouthWest');
set(gca,'FontSize',13);
legend boxoff
xlabel('Threshold of Power','FontSize',15);
ylabel('Area Under Curve','FontSize',15);
xlim([1/limN 1]);
ylim([0 1]);
titleStr = strcat('AUC of Performance Profiles for High-Dimensional Settings');
title(titleStr,'FontSize',14);
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)
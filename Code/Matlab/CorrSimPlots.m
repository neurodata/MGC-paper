function []=CorrSimPlots(total,pre1,pre2)
% Used to plot figure 1-8 used in tex. Run like

% total is usually 20.
% pre1 specifies the location to load data.
% pre2 specifies the location to save pictures.

if nargin<1
    total=20; % Usually 20, but can be changed in case of new simulations
end
if nargin<2
    pre1='../../Data/'; % The folder to locate data
end
if nargin<3
    pre2='../../Figures/Fig'; % The folder to save figures
end

%% Set colors
map1=zeros(7,3);
cm = 3;

switch cm
    case 1
        cmap = brewermap(8,'Dark2');
    case 2
        cmap = brewermap(8,'Set2');
    case 3
        gr = [0,1,0];
        ma = [1,0,1];
        cy = [0,1,1];
        cmap(1,:) = gr;
        cmap(2,:) = ma;
        cmap(3,:) = cy;
        cmap(7,:) = [0 0 0];
    case 4
        cmap(1,:) = [166,206,227]/255;
        cmap(2,:) = [31,120,180]/255;
        cmap(3,:) = [178,223,138]/255;
    case 5
        cmap(1,:) = [102,194,165]/255;
        cmap(2,:) = [ 252,141,98]/255;
        cmap(3,:) = [141,160,203]/255;
end

mcorr = cmap(1,:);
dcorr = cmap(2,:);
mante = cmap(3,:);
HHG   = [0,0,0];

map1(1,:)=mcorr; map1(4,:)=mcorr; % The color for MGC{mcorr} and global mcorr.
map1(2,:)=dcorr; map1(5,:)=dcorr; % The color for MGC{dcorr} and global dcorr.
map1(3,:)=mante; map1(6,:)=mante; % The color for MGC{Mantel} and global Mantel.
map1(7,:)=HHG; % The color for HHG
map2 = brewermap(128,'GnBu'); % brewmap

%figure1-4
figNumber='1';
figure('units','normalized','position',[0 0 1 1])
set(groot,'defaultAxesColorOrder',map1);
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    %
    plot(numRange,power1,'.-',numRange,power2,'.-',numRange,power3,'.-',numRange,power4,'.:',numRange,power5,'.:',numRange,power6,'.:',numRange,power7,'.--','LineWidth',3);
    xlim([numRange(1) numRange(end)]);
    ylim([0 1]);
    if j~=16 % Remove x&y axis ticks except type 16, which is at the left bottom
        set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
    end
    title(titlechar);
end
xlabel('Sample Size','position',[-200 -0.2],'FontSize',20);
ylabel('Empirical Testing Power','position',[-515 3],'FontSize',20);
h=suptitle('Testing Powers of 20 Simulated Dependencies for Dimension 1 with Increasing Sample Size');
set(h,'FontSize',20,'FontWeight','normal');
lgdPosition = [0.03, 0.87, .07, .07]; %Legend Position
h=legend('MGC\{mcorr\}','MGC\{dcorr\}','MGC\{Mantel\}','mcorr','dcorr','Mantel','HHG','Location',lgdPosition);
set(h,'FontSize',12);
%
F.fname=[strcat(pre2, figNumber)]; %, '_', num2str(cm)];
F.wh=[8 4]*2;
print_fig(gcf,F)

%%
figNumber='2';
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    K=n;kmin=1;thres=0.8;
    ind=[find(max(power1,[],1)>=thres,1) lim];
    lim=min(ind);
    kmin=2;
    ph=power1All(kmin:numRange(lim),kmin:numRange(lim),lim)';
    imagesc(ph);
    %set(gca,'YDir','normal')
    colormap(map2)
    caxis([0 thres])
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    title(titlechar);
end
xlabel('Neighborhood Choice of X','position',[-220 120],'FontSize',20);
ylabel('Neighborhood Choice of Y','position',[-540 -200],'FontSize',20);
colorbar
tstring=' of mcorr ';
h=suptitle(strcat('Testing Powers of All Local Tests',tstring, ' for Dimension 1'));
set(h,'FontSize',20,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 4]*2;
print_fig(gcf,F)

%% %performance profile
figNumber='3';
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
    ind=[find(power1>=thres,1) find(power4>=thres,1) find(power5>=thres,1) find(power6>=thres,1) find(power7>=thres,1) lim];
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
plot(xaxis,profile(1,:),'.-',xaxis,profile(2,:),'.-',xaxis,profile(3,:),'.-',xaxis,profile(4,:),'.:',xaxis, profile(5,:),'.:',xaxis,profile(6,:),'.:',xaxis,profile(7,:),'.--','LineWidth',2);
h=legend(strcat('MGC\{mcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{dcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{Mantel\}, AUC=', num2str(sumP(3))),strcat('mcorr, AUC=', num2str(sumP(4))),strcat('dcorr, AUC=', num2str(sumP(5))),strcat('Mantel, AUC=', num2str(sumP(6))),strcat('HHG, AUC=', num2str(sumP(7))),'Location','SouthEast');
set(h,'FontSize',12);
xlabel('Difference with the Best Method','FontSize',16);
ylabel('Relative Performance','FontSize',16);
ylim([0 1]);
titleStr = strcat('Performance Profiles for Dimension 1');
title(titleStr,'FontSize',12);
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%%performance profile AUC for n
figNumber='4';
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
        ind=[find(power1>=thres,1) find(power4>=thres,1) find(power5>=thres,1) find(power6>=thres,1) find(power7>=thres,1) lim];
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
plot(xaxis,sumP(1,:),'.-',xaxis,sumP(2,:),'.-',xaxis,sumP(3,:),'.-',xaxis, sumP(4,:),'.:',xaxis,sumP(5,:),'.:',xaxis, sumP(6,:),'.:',xaxis,sumP(7,:),'.--','LineWidth',2);
h=legend('MGC\{mcorr\}','MGC\{dcorr\}','MGC\{Mantel\}','mcorr','dcorr','Mantel','HHG','Location','SouthEast');
set(h,'FontSize',12);
xlabel('Threshold of Power','FontSize',16);
ylabel('Area Under Curve','FontSize',16);
ylim([0 1]);
titleStr = strcat('Area Under Curve of Performance Profiles for Dimension 1');
title(titleStr,'FontSize',12);
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)


%Plot 5-8
figNumber='5';
figure('units','normalized','position',[0 0 1 1])
set(groot,'defaultAxesColorOrder',map1)
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    numRange=dimRange;
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    plot(numRange,power1,'.-',numRange,power2,'.-',numRange,power3,'.-',numRange,power4,'.:',numRange,power5,'.:',numRange,power6,'.:',numRange,power7,'.--','LineWidth',2);
    xlim([numRange(1) numRange(end)]);
    ylim([0 1]);
    if j~=16 % Remove x&y axis ticks except type 16, which is at the left bottom
        %set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
    end
    title(titlechar);
end
xlabel('Dimension','position',[-200 -0.2],'FontSize',20);
ylabel('Empirical Testing Power','position',[-515 3],'FontSize',20);
h=suptitle('Testing Powers of 20 Simulated Dependencies for Increasing Dimension with Fixed Sample Size');
set(h,'FontSize',20,'FontWeight','normal');
lgdPosition = [0.03, 0.87, .07, .07]; %Legend Position
h=legend('MGC\{mcorr\}','MGC\{dcorr\}','MGC\{Mantel\}','mcorr','dcorr','Mantel','HHG','Location',lgdPosition);
set(h,'FontSize',12);
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 4]*2;
print_fig(gcf,F)


figNumber='6';
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    kmin=1;thres=0.5;
    ind=[find(max(power1,[],1)>=thres,1,'last'),1];
    lim=max(ind);
    ph=power1All(kmin:n,kmin:n,lim)';
%     if max(max(ph))>thres
%         ph=ph/max(max(ph))*thres; % in cas
%     end
    imagesc(ph);
    %set(gca,'YDir','normal')
    colormap(map2)
    caxis([0 thres])
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    title(titlechar);
end
xlabel('Neighborhood Choice of X','position',[-220 120],'FontSize',20);
ylabel('Neighborhood Choice of Y','position',[-540 -200],'FontSize',20);
colorbar
tstring=' of mcorr ';
h=suptitle(strcat('Testing Powers of All Local Tests',tstring,' for Increasing Dimension'));
set(h,'FontSize',20,'FontWeight','normal');
%
F.fname=strcat(pre2, figNumber);
F.wh=[8 4]*2;
print_fig(gcf,F)

%%%performance profile
figNumber='7';
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
    ind=[find(power1>=thres,1,'last') find(power4>=thres,1,'last') find(power5>=thres,1,'last') find(power6>=thres,1,'last') find(power7>=thres,1,'last') 1];
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
plot(xaxis,profile(1,:),'.-',xaxis,profile(2,:),'.-',xaxis,profile(3,:),'.-',xaxis,profile(4,:),'.:',xaxis, profile(5,:),'.:',xaxis,profile(6,:),'.:',xaxis,profile(7,:),'.--','LineWidth',2);
h=legend(strcat('MGC\{mcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{dcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{Mantel\}, AUC=', num2str(sumP(3))),strcat('mcorr, AUC=', num2str(sumP(4))),strcat('dcorr, AUC=', num2str(sumP(5))),strcat('Mantel, AUC=', num2str(sumP(6))),strcat('HHG, AUC=', num2str(sumP(7))),'Location','SouthEast');
set(h,'FontSize',12);
xlabel('Difference with the Best Method','FontSize',16);
ylabel('Relative Performance','FontSize',16);
ylim([0 1]);
titleStr = strcat('Performance Profiles for Increasing Dimension');
title(titleStr,'FontSize',12);
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%%performance profile
figNumber='8';
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
        ind=[find(power1>=thres,1,'last') find(power4>=thres,1,'last') find(power5>=thres,1,'last') find(power6>=thres,1,'last') find(power7>=thres,1,'last') 1];
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
plot(xaxis,sumP(1,:),'.-',xaxis,sumP(2,:),'.-',xaxis,sumP(3,:),'.-',xaxis, sumP(4,:),'.:',xaxis,sumP(5,:),'.:',xaxis, sumP(6,:),'.:',xaxis,sumP(7,:),'.--','LineWidth',2);
h=legend('MGC\{mcorr\}','MGC\{dcorr\}','MGC\{Mantel\}','mcorr','dcorr','Mantel','HHG','Location','SouthEast');
set(h,'FontSize',12);
xlabel('Threshold of Power','FontSize',16);
ylabel('Area Under Curve','FontSize',16);
ylim([0 1]);
titleStr = strcat('Area Under Curve of Performance Profiles for Increasing Dimension');
title(titleStr,'FontSize',12);
%
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)
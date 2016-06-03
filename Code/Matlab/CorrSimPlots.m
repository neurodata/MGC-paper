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
    pre2='../../Draft/Figures/Fig'; % The folder to save figures
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

dcorr = cmap(1,:);
mcorr = cmap(2,:);
mante = cmap(3,:);
HHG   = [0,0,0];

map1(1,:)=dcorr; map1(4,:)=dcorr; % The color for MGC{dcorr} and global dcorr.
map1(2,:)=mcorr; map1(5,:)=mcorr; % The color for MGC{mcorr} and global mcorr.
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
    set(gca,'FontSize',14);
    title(titlechar,'FontSize',14);
end
xlabel('Sample Size','position',[-200 -0.2],'FontSize',20);
ylabel('Empirical Testing Power','position',[-515 3],'FontSize',20);
h=suptitle('Testing Powers for 20 Simulated 1-Dimensional Settings');
set(h,'FontSize',24,'FontWeight','normal');
lgdPosition = [0.05, 0.87, .05, .05]; %Legend Position
h=legend('MGC\{dcorr\}','MGC\{mcorr\}','MGC\{Mantel\}','dcorr','mcorr','Mantel','HHG','Location',lgdPosition);
legend boxoff
set(h,'FontSize',14);
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
    K=n;kmin=2;thres=0.8;
    ind=[find(max(power2,[],1)>=thres,1) lim];
    lim=min(ind);
    ph=power2All(kmin:numRange(lim),kmin:numRange(lim),lim)';
 tt=find(sum(ph,2)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(tt:end,:)=repmat(ph(tt-1,:),numRange(lim)-tt,1);
    end
    tt=find(sum(ph,1)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(:,tt:end)=repmat(ph(:,tt-1),1,numRange(lim)-tt);
    end
    imagesc(ph);
    set(gca,'YDir','normal')
    colormap(map2)
    caxis([0 thres])
    set(gca,'FontSize',12);
%     set(gca,'XTick',[]); % Remove x axis ticks
%     set(gca,'YTick',[]); % Remove y axis ticks
    title(titlechar,'FontSize',14);
end
xlabel('Number of Neighbors for X','position',[-210 -20],'FontSize',20);
ylabel('Number of Neighbors for Y','position',[-540 300],'FontSize',20);
colorbar
tstring=' of mcorr ';
h=suptitle(strcat('Testing Powers of All Local Correlations for 1-Dimensional Simulations'));
set(h,'FontSize',24,'FontWeight','normal');
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
plot(xaxis,profile(1,:),'.-',xaxis,profile(2,:),'.-',xaxis,profile(3,:),'.-',xaxis,profile(4,:),'.:',xaxis, profile(5,:),'.:',xaxis,profile(6,:),'.:',xaxis,profile(7,:),'.--','LineWidth',2);
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
set(gca,'FontSize',13);
legend(strcat('MGC\{dcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{mcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{Mantel\}, AUC=', num2str(sumP(3))),strcat('dcorr, AUC=', num2str(sumP(4))),strcat('mcorr, AUC=', num2str(sumP(5))),strcat('Mantel, AUC=', num2str(sumP(6))),strcat('HHG, AUC=', num2str(sumP(7))),'Location','SouthEast');
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
plot(xaxis,sumP(1,:),'.-',xaxis,sumP(2,:),'.-',xaxis,sumP(3,:),'.-',xaxis, sumP(4,:),'.:',xaxis,sumP(5,:),'.:',xaxis, sumP(6,:),'.:',xaxis,sumP(7,:),'.--','LineWidth',2);
set(gca,'XTick',[]); % Remove x axis ticks
set(gca,'YTick',[]); % Remove y axis ticks
legend('MGC\{dcorr\}','MGC\{mcorr\}','MGC\{Mantel\}','dcorr','mcorr','Mantel','HHG','Location','SouthWest');
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
    set(gca,'FontSize',12);
    title(titlechar,'FontSize',14);
end
xlabel('Dimension','position',[-200 -0.2],'FontSize',20);
ylabel('Empirical Testing Power','position',[-515 3],'FontSize',20);
h=suptitle('Testing Powers for 20 Simulated High-Dimensional Settings');
set(h,'FontSize',24,'FontWeight','normal');
lgdPosition = [0.05, 0.87, .05, .05]; %Legend Position
h=legend('MGC\{dcorr\}','MGC\{mcorr\}','MGC\{Mantel\}','dcorr','mcorr','Mantel','HHG','Location',lgdPosition);
set(h,'FontSize',14);
legend boxoff
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
    kmin=2;thres=0.5;
    ind=[find(max(power2,[],1)>=thres,1,'last'),1];
    lim=max(ind);
    ph=power2All(kmin:n,kmin:n,lim)';
 tt=find(sum(ph,2)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(tt:end,:)=repmat(ph(tt-1,:),n-tt,1);
    end
    tt=find(sum(ph,1)==0,1,'first');
    if isempty(tt)==false && tt~=1;
        ph(:,tt:end)=repmat(ph(:,tt-1),1,n-tt);
    end
%     if max(max(ph))>thres
%         ph=ph/max(max(ph))*thres; % in cas
%     end
    imagesc(ph);
    set(gca,'YDir','normal')
    colormap(map2)
    caxis([0 thres])
    if j~=16
    set(gca,'XTick',[]); % Remove x axis ticks
    set(gca,'YTick',[]); % Remove y axis ticks
    end
    set(gca,'FontSize',14);
    title(titlechar,'FontSize',14);
end
xlabel('Number of Neighbors for X','position',[-210 -20],'FontSize',20);
ylabel('Number of Neighbors for Y','position',[-540 300],'FontSize',20);
colorbar
h=suptitle(strcat('Testing Powers of All Local Correlations for High-Dimensional Simulations'));
set(h,'FontSize',24,'FontWeight','normal');
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
plot(xaxis,profile(1,:),'.-',xaxis,profile(2,:),'.-',xaxis,profile(3,:),'.-',xaxis,profile(4,:),'.:',xaxis, profile(5,:),'.:',xaxis,profile(6,:),'.:',xaxis,profile(7,:),'.--','LineWidth',2);
legend(strcat('MGC\{dcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{mcorr\}, AUC=', num2str(sumP(1))),strcat('MGC\{Mantel\}, AUC=', num2str(sumP(3))),strcat('dcorr, AUC=', num2str(sumP(4))),strcat('mcorr, AUC=', num2str(sumP(5))),strcat('Mantel, AUC=', num2str(sumP(6))),strcat('HHG, AUC=', num2str(sumP(7))),'Location','SouthEast');
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
plot(xaxis,sumP(1,:),'.-',xaxis,sumP(2,:),'.-',xaxis,sumP(3,:),'.-',xaxis, sumP(4,:),'.:',xaxis,sumP(5,:),'.:',xaxis, sumP(6,:),'.:',xaxis,sumP(7,:),'.--','LineWidth',2);
legend('MGC\{dcorr\}','MGC\{mcorr\}','MGC\{Mantel\}','dcorr','mcorr','Mantel','HHG','Location','SouthWest');
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


%
map1(1,:)=dcorr;
map1(2,:)=mcorr; 
map1(3,:)=mcorr; 
map1(4,:)=HHG; 
set(groot,'defaultAxesColorOrder',map1);

figNumber='9';
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
plot(x,p1(1,:),'.-',x,p1(3,:),'.-',x,p1(2,:),'.:',x,p1(4,:),'.--','LineWidth',2);
legend('Estimated MGC', 'True MGC','Global mcorr','HHG','Location','SouthWest');
set(gca,'FontSize',14);
legend boxoff
xlabel('Function Type','FontSize',15);
ylabel('Testing Power','FontSize',15);
ylim([0,1]);
title('1-Dimensional Simulations at n=60','FontSize',18);
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

figNumber='10';
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
plot(x,p1(1,:),'.-',x,p1(3,:),'.-',x,p1(2,:),'.:',x,p1(4,:),'.--','LineWidth',2);
%h=legend('Estimated MGC', 'True MGC','Global Mcorr','Location','SouthWest');
%set(h,'FontSize',12);
set(gca,'FontSize',14);
xlabel('Function Type','FontSize',15);
ylabel('Testing Power','FontSize',15);
ylim([0,1]);
title('High-Dimensional Simulations at n=100','FontSize',18);
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%% Outlier Model
seq=[0.3,0.5,0.7];
for i=1:length(seq)
    figNumber=strcat('Out',num2str(i));
    filename=strcat(pre1,'CorrIndTestType0N100Dim1P',num2str(seq(i)),'.mat');
    load(filename)
    figure
    imagesc(power2All(2:end,2:end));
    set(gca,'YDir','normal')
    set(gca,'FontSize',14);
    colormap(map2)
    caxis([0.5 1])
    if i==2;
        colorbar();
    else
        set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
    end
    title(strcat('Local Testing Power at p=',num2str(seq(i))),'FontSize',15);
    F.fname=strcat(pre2, figNumber);
    F.wh=[3 2.5]*2;
    print_fig(gcf,F)
end
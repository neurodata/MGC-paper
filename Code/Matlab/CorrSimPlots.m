function []=CorrSimPlots(optionA, total,pre1,pre2)
% Used to plot figure 1-8 used in tex. Run like
% optionA=1;
% CorrSimPlots(optionA)
% Note that there are still some problems for fig2&6 plot size, so these
% two should be manually saved...


% optionA can be 1 or 2. Use 1 for paper figures, which include MGC by mcorr and all global test; 
% use 2 for appendix figures, which will further include MGC by dcorr/Mantel.
%
% total is usually 20.
% pre1 specifies the location to load data.
% pre2 specifies the location to save pictures.

if nargin<1
    optionA=1; % Use 1 for paper figures, which include MGC by mcorr and all global test; use 2 for appendix figures, which will further include MGC by dcorr/Mantel.
end
if nargin<2
    total=20;
end
if nargin<3
    pre1='../../Data/'; % The folder to locate data
    %pre1='News_1/';
end
if nargin<4
    pre2='../../Figures/JovoFig'; % The folder to save figures
    %pre2='News_1/Fig';
end

%figure1-4
figNumber='1';
if optionA~=1
    figNumber=strcat(figNumber,'b');
end
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);power1L=zeros(length(numRange),1);power2L=zeros(length(numRange),1);power3L=zeros(length(numRange),1);power1M=zeros(length(numRange),1);power2M=zeros(length(numRange),1);power3M=zeros(length(numRange),1);
    for i=1:length(numRange)
        power1L(i)=power1(numRange(i),numRange(i),i);
        power2L(i)=power2(numRange(i),numRange(i),i);
        power3L(i)=power3(numRange(i),numRange(i),i);
        power1M(i)=max(max(power1(2:end,2:end,i)));
        power2M(i)=max(max(power2(2:end,2:end,i)));
        power3M(i)=max(max(power3(2:end,2:end,i)));
    end
    switch optionA
        case 1
            plot(numRange,power1M,'ro-',numRange,power1L,'r.: ',numRange,power2L,'b.:',numRange,power3L,'c.:',numRange,power4,'g.:','LineWidth',2);           
        case 2
            plot(numRange,power1M,'ro-',numRange,power2M,'bx-',numRange,power3M,'c+-',numRange,power1L,'r.:',numRange,power2L,'b.:',numRange,power3L,'c.:',numRange,power4,'g.:','LineWidth',2);
    end
    xlim([numRange(1) numRange(end)]);
    ylim([0 1]);
    title(titlechar);
end
xlabel('Sample Size','position',[-200 -0.2],'FontSize',20);
ylabel('Empirical Testing Power','position',[-520 3],'FontSize',20);
h=suptitle('Testing Powers of 20 Simulated Dependencies for Dimension 1 with Increasing Sample Size');
set(h,'FontSize',20,'FontWeight','normal');
lgdPosition = [0.03, 0.87, .07, .07]; %Legend Position
switch optionA
    case 1
        h=legend('MGC','mcorr','dcorr','Mantel','HHG','Location',lgdPosition);    
    case 2
        h=legend('MGC by mcorr','MGC by dcorr','MGC by Mantel','mcorr','dcorr','Mantel','HHG','Location',lgdPosition);
end
set(h,'FontSize',12);
%
F.fname=[strcat(pre2, figNumber)];
F.wh=[8 4]*2;
print_fig(gcf,F)

%
figNumber='2';
if optionA~=1
    figNumber=strcat(figNumber,'b');
end
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    if optionA==1
    else
        power1=power3;
    end
    for i=1:length(numRange)
        p1(i,:)=max(power1(2:end,:,i),[],1);
        %         p2(i,:)=max(power2(2:end,:,i),[],1);
        %         %         p1(i,:)=power1(i,numRange(i),:);%
        %         %         p2(i,:)=power2(i,numRange(i),:);%
        %         %         p3(i,:)=power3(i,numRange(i),:);%check
        %         power1L(i)=power1(numRange(i),numRange(i),i);
        %         power2L(i)=power2(numRange(i),numRange(i),i);
    end
    K=n;kmin=1;thres=0.8;
    %ind=[find(max(power1,[],2)>thres,1) find(max(power2,[],2)>thres,1) find(power3>thres,1) find(power4>thres,1) lim];
    ind=[find(max(p1,[],2)>=thres,1) lim];
    lim=min(ind);
    if numRange(lim)>50
        c=2;
        lim=floor(lim/2);
    else
        c=1;
        kmin=2;
    end
    xaxis=kmin:numRange(lim);
    yaxis=kmin:numRange(lim);
    [X,Y]=meshgrid(c*xaxis,c*yaxis);
    ph=power1(c*xaxis,c*yaxis,c*lim)';
    surf(X,Y,ph);
    view(2)
    caxis([0 thres])
    xlim([c*kmin c*numRange(lim)]);
    ylim([c*kmin c*numRange(lim)]);
    title(titlechar);
end
xlabel('Neighborhood Choice of X','position',[-200 -25],'FontSize',20);
ylabel('Neighborhood Choice of Y','position',[-530 300],'FontSize',20);
colorbar
if optionA==1
    tstring=' by mcorr ';
else
    tstring=' by Mantel ';
end
h=suptitle(strcat('Testing Power of Multiscale Graph Dependency',tstring, ' for Dimension 1'));
set(h,'FontSize',20,'FontWeight','normal');
%
% F.fname=[strcat(pre2, figNumber)];
% F.wh=[8 4]*2;
% print_fig(gcf,F)

%%%performance profile
figNumber='3';
if optionA~=1
    figNumber=strcat(figNumber,'b');
end
figure
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
profile=zeros(7,length(xaxis));
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    for i=1:length(numRange)
        power1L(i)=power1(numRange(i),numRange(i),i);
        power2L(i)=power2(numRange(i),numRange(i),i);
        power3L(i)=power3(numRange(i),numRange(i),i);
        power1M(i)=max(max(power1(2:end,2:end,i)));
        power2M(i)=max(max(power2(2:end,2:end,i)));
        power3M(i)=max(max(power3(2:end,2:end,i)));
    end
    thres=0.8;
    ind=[find(power1M>=thres,1) find(power2L>=thres,1) find(power3L>=thres,1) find(power4>=thres,1) lim];
    pos=min(ind);
    switch optionA
        case 1
            power=[power1M(pos), power1L(pos), power2L(pos), power3L(pos),power4(pos)];            
        case 2
            power=[power1M(pos), power1L(pos), power2L(pos), power3L(pos),power4(pos),power2M(pos),power3M(pos)];
    end
    pmax=max(power);
    for k=1:length(power)
        tmp=(pmax-power(k));
        tmpInd=ceil(tmp/interval)+1;
        profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
    end
end
profile=profile./total;
sumP=ceil(mean(profile,2)*1000)/1000;
switch optionA
    case 1
        plot(xaxis,profile(1,:),'r.-',xaxis,profile(2,:),'r.:',xaxis, profile(3,:),'b.:',xaxis,profile(4,:),'c.:',xaxis,profile(5,:),'g.:','LineWidth',2);
        h=legend(strcat('MGC, AUC=', num2str(sumP(1))),strcat('mcorr, AUC=', num2str(sumP(2))),strcat('dcorr, AUC=', num2str(sumP(3))),strcat('Mantel, AUC=', num2str(sumP(4))),strcat('HHG, AUC=', num2str(sumP(5))),'Location','SouthEast');
    case 2
        plot(xaxis,profile(1,:),'r.-',xaxis,profile(6,:),'b.-',xaxis,profile(7,:),'c.-',xaxis,profile(2,:),'r.:',xaxis, profile(3,:),'b.:',xaxis,profile(4,:),'c.:',xaxis,profile(5,:),'g.:','LineWidth',2);
        h=legend(strcat('MGC by mcorr, AUC=', num2str(sumP(1))),strcat('MGC by dcorr, AUC=', num2str(sumP(6))),strcat('MGC by Mantel, AUC=', num2str(sumP(7))),strcat('mcorr, AUC=', num2str(sumP(2))),strcat('dcorr, AUC=', num2str(sumP(3))),strcat('Mantel, AUC=', num2str(sumP(4))),strcat('HHG, AUC=', num2str(sumP(5))),'Location','SouthEast');
end
set(h,'FontSize',12);
xlabel('Difference with the Best Method','FontSize',13);
ylabel('Relative Performance','FontSize',13);
ylim([0 1]);
titleStr = strcat('Performance Profiles for Dimension 1');
title(titleStr,'FontSize',12);
%
F.fname=[strcat(pre2, figNumber)];
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%%performance profile AUC for n
figNumber='4';
if optionA~=1
    figNumber=strcat(figNumber,'b');
end
figure
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
        for i=1:length(numRange)
            power1L(i)=power1(numRange(i),numRange(i),i);
            power2L(i)=power2(numRange(i),numRange(i),i);
            power3L(i)=power3(numRange(i),numRange(i),i);
            power1M(i)=max(max(power1(2:end,2:end,i)));
            power2M(i)=max(max(power2(2:end,2:end,i)));
            power3M(i)=max(max(power3(2:end,2:end,i)));
        end
        ind=[find(power1M>=thres,1) find(power2L>=thres,1) find(power3L>=thres,1) find(power4>=thres,1) lim];
        pos=min(ind);
        switch optionA
            case 1
                power=[power1M(pos), power1L(pos), power2L(pos), power3L(pos),power4(pos)];               
            case 2
                power=[power1M(pos), power1L(pos), power2L(pos), power3L(pos),power4(pos),power2M(pos),power3M(pos)];
        end
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
switch optionA
    case 1
        plot(xaxis,sumP(1,:),'r.-',xaxis, sumP(2,:),'r.:',xaxis,sumP(3,:),'b.:',xaxis, sumP(4,:),'c.:',xaxis,sumP(5,:),'g.:','LineWidth',2);
        h=legend('MGC','mcorr','dcorr','Mantel','HHG','Location','SouthEast');
    case 2
        plot(xaxis,sumP(1,:),'r.-',xaxis,sumP(6,:),'b.-',xaxis,sumP(7,:),'c.-',xaxis, sumP(2,:),'r.:',xaxis,sumP(3,:),'b.:',xaxis, sumP(4,:),'c.:',xaxis,sumP(5,:),'g.:','LineWidth',2);
        h=legend('MGC by mcorr','MGC by dcorr','MGC by Mantel','mcorr','dcorr','Mantel','HHG','Location','SouthEast');
end
set(h,'FontSize',12);
xlabel('Threshold of Power','FontSize',16);
ylabel('Area Under Curve','FontSize',16);
ylim([0 1]);
titleStr = strcat('Area Under Curve of Performance Profiles for Dimension 1');
title(titleStr,'FontSize',12);
%
F.fname=[strcat(pre2, figNumber)];
F.wh=[3 2.5]*2;
print_fig(gcf,F)


%Plot 5-8
figNumber='5';
if optionA~=1
    figNumber=strcat(figNumber,'b');
end
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    numRange=dimRange;
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    power1M=zeros(length(dimRange),1);power2M=zeros(length(dimRange),1);power3M=zeros(length(dimRange),1);power1L=zeros(length(dimRange),1);power2L=zeros(length(dimRange),1);power3L=zeros(length(dimRange),1);
    for i=1:length(dimRange)
        power1L(i)=power1(n,n,i);
        power2L(i)=power2(n,n,i);
        power3L(i)=power3(n,n,i);
        power1M(i)=max(max(power1(2:end,2:end,i)));
        power2M(i)=max(max(power2(2:end,2:end,i)));      
        power3M(i)=max(max(power3(2:end,2:end,i)));
    end
    switch optionA
        case 1
            plot(numRange,power1M,'ro-',numRange,power1L,'r.:',numRange,power2L,'b.:',numRange,power3L,'c.:',numRange,power4,'g.:','LineWidth',2);           
        case 2
            plot(numRange,power1M,'ro-',numRange,power2M,'bx-',numRange,power3M,'c+-',numRange,power1L,'r.:',numRange,power2L,'b.:',numRange,power3L,'c.:',numRange,power4,'g.:','LineWidth',2);
    end
    xlim([numRange(1) numRange(end)]);
    ylim([0 1]);
    title(titlechar);
end
xlabel('Dimension','position',[-200 -0.2],'FontSize',20);
ylabel('Empirical Testing Power','position',[-520 3],'FontSize',20);
h=suptitle('Testing Powers of 20 Simulated Dependencies for Increasing Dimension with Fixed Sample Size');
set(h,'FontSize',20,'FontWeight','normal');
lgdPosition = [0.03, 0.87, .07, .07]; %Legend Position
switch optionA
    case 1
        h=legend('MGC','mcorr','dcorr','Mantel','HHG','Location',lgdPosition);
    case 2
        h=legend('MGC by mcorr','MGC by dcorr','MGC by Mantel','mcorr','dcorr','Mantel','HHG','Location',lgdPosition);
end
set(h,'FontSize',12);
%
F.fname=[strcat(pre2, figNumber)];
F.wh=[8 4]*2;
print_fig(gcf,F)

figNumber='6';
if optionA~=1
    figNumber=strcat(figNumber,'b');
end
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    p1=zeros(length(dimRange),n);
    if optionA==1;
    else
        power1=power3;
    end
    for i=1:length(dimRange)
        p1(i,:)=max(power1(2:end,:,i),[],1);
        %         p2(i,:)=max(power2(2:end,:,i),[],1);
        %         power1L(i)=power1(end,end,i);
        %         power2L(i)=power2(end,end,i);
    end
    K=n;kmin=1;thres=0.5;
    ind=[find(max(p1,[],2)>=thres,1,'last'),1];
    lim=max(ind);
    numLim=50;
    xaxis=kmin:numLim;
    yaxis=kmin:numLim;
    [X,Y]=meshgrid(2*xaxis,2*yaxis);
    ph=power1(2*xaxis,2*yaxis,lim)';
    surf(X,Y,ph);
    view(2)
    caxis([0 thres])
    xlim([2*kmin 2*numLim]);
    ylim([2*kmin 2*numLim]);
    title(titlechar);
end
xlabel('Neighborhood Choice of X','position',[-200 -25],'FontSize',20);
ylabel('Neighborhood Choice of Y','position',[-530 300],'FontSize',20);
colorbar
if optionA==1
    tstring=' by mcorr ';
else
    tstring=' by Mantel ';
end
h=suptitle(strcat('Testing Power of Multiscale Graph Dependency',tstring,' for Increasing Dimension'));
set(h,'FontSize',20,'FontWeight','normal');
%
% F.fname=[strcat(pre2, figNumber)];
% F.wh=[8 4]*2;
% print_fig(gcf,F)

%%%performance profile
figNumber='7';
if optionA~=1
    figNumber=strcat(figNumber,'b');
end
figure
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
profile=zeros(7,length(xaxis));
%load data
for j=1:total
    filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    power1M=zeros(length(dimRange),1);power2M=zeros(length(dimRange),1);power3M=zeros(length(dimRange),1);power1L=zeros(length(dimRange),1);power2L=zeros(length(dimRange),1);power3L=zeros(length(dimRange),1);
    for i=1:length(dimRange)
        power1L(i)=power1(n,n,i);
        power2L(i)=power2(n,n,i);
        power3L(i)=power3(n,n,i);
        power1M(i)=max(max(power1(2:end,2:end,i)));
        power2M(i)=max(max(power2(2:end,2:end,i)));
        power3M(i)=max(max(power3(2:end,2:end,i)));
    end
    thres=0.5;
    ind=[find(power1M>=thres,1,'last') find(power2L>=thres,1,'last') find(power3L>=thres,1,'last') find(power4>=thres,1,'last') 1];
    pos=max(ind);
    switch optionA
        case 1
            power=[power1M(pos), power1L(pos), power2L(pos), power3L(pos),power4(pos)];
        case 2
            power=[power1M(pos), power1L(pos), power2L(pos), power3L(pos),power4(pos),power2M(pos),power3M(pos)];
    end
    pmax=max(power);
    for k=1:length(power)
        tmp=(pmax-power(k));
        tmpInd=ceil(tmp/interval)+1;
        profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
    end
end
profile=profile./total;
sumP=ceil(mean(profile,2)*1000)/1000;
switch optionA
    case 1
        plot(xaxis,profile(1,:),'r.-',xaxis,profile(2,:),'r.:',xaxis, profile(3,:),'b.:',xaxis,profile(4,:),'c.:',xaxis,profile(5,:),'g.:','LineWidth',2);
        h=legend(strcat('MGC, AUC=', num2str(sumP(1))),strcat('mcorr, AUC=', num2str(sumP(2))),strcat('dcorr, AUC=', num2str(sumP(3))),strcat('Mantel, AUC=', num2str(sumP(4))),strcat('HHG, AUC=', num2str(sumP(5))),'Location','SouthEast');
    case 2
        plot(xaxis,profile(1,:),'r.-',xaxis,profile(6,:),'b.-',xaxis,profile(7,:),'c.-',xaxis,profile(2,:),'r.:',xaxis, profile(3,:),'b.:',xaxis,profile(4,:),'c.:',xaxis,profile(5,:),'g.:','LineWidth',2);
        h=legend(strcat('MGC by mcorr, AUC=', num2str(sumP(1))),strcat('MGC by dcorr, AUC=', num2str(sumP(6))),strcat('MGC by Mantel, AUC=', num2str(sumP(7))),strcat('mcorr AUC=', num2str(sumP(2))),strcat('dcorr, AUC=', num2str(sumP(3))),strcat('Mantel, AUC=', num2str(sumP(4))),strcat('HHG, AUC=', num2str(sumP(5))),'Location','SouthEast');
end
set(h,'FontSize',12);
xlabel('Difference with the Best Method','FontSize',16);
ylabel('Relative Performance','FontSize',16);
ylim([0 1]);
titleStr = strcat('Performance Profiles for Increasing Dimension');
title(titleStr,'FontSize',12);
%
F.fname=[strcat(pre2, figNumber)];
F.wh=[3 2.5]*2;
print_fig(gcf,F)

%%%performance profile
figNumber='8';
if optionA~=1
    figNumber=strcat(figNumber,'b');
end
figure
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
        power1M=zeros(length(dimRange),1);power2M=zeros(length(dimRange),1);power3M=zeros(length(dimRange),1);power1L=zeros(length(dimRange),1);power2L=zeros(length(dimRange),1);power3L=zeros(length(dimRange),1);
        for i=1:length(dimRange)
            power1L(i)=power1(n,n,i);
            power2L(i)=power2(n,n,i);
            power3L(i)=power3(n,n,i);
            power1M(i)=max(max(power1(2:end,2:end,i)));
            power2M(i)=max(max(power2(2:end,2:end,i)));
            power3M(i)=max(max(power3(2:end,2:end,i)));
        end
        ind=[find(power1M>=thres,1,'last') find(power2L>=thres,1,'last') find(power3L>=thres,1,'last') find(power4>=thres,1,'last') 1];
        pos=max(ind);
        switch optionA
            case 1
                power=[power1M(pos), power1L(pos), power2L(pos), power3L(pos),power4(pos)];            
            case 2
                power=[power1M(pos), power1L(pos), power2L(pos), power3L(pos),power4(pos),power2M(pos),power3M(pos)];
        end
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
switch optionA
    case 1
        plot(xaxis,sumP(1,:),'r.-',xaxis, sumP(2,:),'r.:',xaxis,sumP(3,:),'b.:',xaxis, sumP(4,:),'c.:',xaxis,sumP(5,:),'g.:','LineWidth',2);
        h=legend('MGC','mcorr','dcorr','Mantel','HHG','Location','SouthEast');
    case 2
        plot(xaxis,sumP(1,:),'r.-',xaxis,sumP(6,:),'b.-',xaxis,sumP(7,:),'c.-',xaxis, sumP(2,:),'r.:',xaxis,sumP(3,:),'b.:',xaxis, sumP(4,:),'c.:',xaxis,sumP(5,:),'g.:','LineWidth',2);
        h=legend('MGC by mcorr','MGC by dcorr','MGC by Mantel','mcorr','dcorr','Mantel','HHG','Location','SouthEast');
end
set(h,'FontSize',12);
xlabel('Threshold of Power','FontSize',16);
ylabel('Area Under Curve','FontSize',16);
ylim([0 1]);
titleStr = strcat('Area Under Curve of Performance Profiles for Increasing Dimension');
title(titleStr,'FontSize',12);
%
F.fname=[strcat(pre2, figNumber)];
F.wh=[3 2.5]*2;
print_fig(gcf,F)
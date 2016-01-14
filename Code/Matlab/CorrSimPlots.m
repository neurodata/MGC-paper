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
    %pre1='Results/';
end
if nargin<4
    pre2='../../Figures/Fig'; % The folder to save figures
    %pre2='Results/Fig';
end

%figure1-4
figNumber='1';
if optionA~=1
    figNumber=strcat(figNumber,'b');
end
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
area1=zeros(total,1);area2=zeros(total,1);area3=zeros(total,1);
area4=zeros(total,1);area5=zeros(total,1);area6=zeros(total,1);
area7=zeros(total,1);
for j=1:total
    filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=CorrSimuTitle(j);
    power1W=zeros(lim,1);power2W=zeros(lim,1);power3W=zeros(lim,1);power7=power4;power4=zeros(lim,1);power5=zeros(lim,1);power6=zeros(lim,1);
    for i=1:lim
        nn=numRange(i);
        power1W(i)=max(max(power1(1:nn,1:nn,i)));
        power2W(i)=max(max(power2(1:nn,1:nn,i)));
        power3W(i)=max(max(power3(1:nn,1:nn,i)));
        power4(i)=power1(nn,nn,i);
        power5(i)=power2(nn,nn,i);
        power6(i)=power3(nn,nn,i);
    end
    power1=power1W;power2=power2W;power3=power3W;
    switch optionA
        case 1
            plot(numRange,power1,'r-',numRange,power4,'r.: ',numRange,power5,'b.:',numRange,power6,'c.:',numRange,power7,'g.:','LineWidth',2);
        case 2
            plot(numRange,power1,'r-',numRange,power2,'bx-',numRange,power3,'c+-',numRange,power4,'r.:',numRange,power5,'b.:',numRange,power6,'c.:',numRange,power7,'g.:','LineWidth',2);
    end
    xlim([numRange(1) numRange(end)]);
    ylim([0 1]);
    title(titlechar);
    area1(j)=mean(power1);
    area2(j)=mean(power2);
    area3(j)=mean(power3);
    area4(j)=mean(power4);
    area5(j)=mean(power5);
    area6(j)=mean(power6);
    area7(j)=mean(power7);
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
% figure

%


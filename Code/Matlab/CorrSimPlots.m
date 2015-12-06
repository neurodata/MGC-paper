function []=CorrSimPlots(option, total)

% Used to plot figure 1-8 used in tex, should be in same folder as all data

clear
option=1;
total=20;
%%figure1-4
figure
s=5;
t=4;
for j=1:total
    filename=strcat('CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=' Data';
    switch type
        case 1
            titlechar=' 1. Linear';
        case 2
            titlechar=' 2. Quadratic';
        case 3
            titlechar=' 3. Cubic';
        case 4
            titlechar=' 4. Sine Period 1/2';
        case 5
            titlechar=' 5. Sine Period 1/8';
        case 6
            titlechar=' 6. X\^(1/4)';
        case 7
            titlechar=' 7. Circle';
        case 8
            titlechar=' 8. Step Function';
        case 9
            titlechar=' 9. Exp(X)';
        case 10
            titlechar=' 10. Uncorrelated Binomial';
        case 11
            titlechar=' 11. W';
        case 12
            titlechar=' 12. Square';
        case 13
            titlechar=' 13. Diamond';
        case 14
            titlechar=' 14. Parabola';
        case 15
            titlechar=' 15. Two Parabolas';
        case 16
            titlechar=' 16. Circle 2';
        case 17
            titlechar=' 17. Independent Clouds';
        case 18
            titlechar=' 18. Joint Normal';
        case 19
            titlechar=' 19. Log(X^2)';
        case 20
            titlechar=' 20. Multiplicative Noise';
    end
    for i=1:length(numRange)
        power1L(i)=power1(numRange(i),numRange(i),i);
        power2L(i)=power2(numRange(i),numRange(i),i);
        power1M(i)=max(max(power1(:,:,i)));
        power2M(i)=max(max(power2(:,:,i)));
    end
    switch option
        case 0
            plot(numRange,power2M,'ro-',numRange,power1M,'bx-',numRange,power2L,'r.:',numRange,power1L,'b.:',numRange,power3,'g.-','LineWidth',1);
        case 1
            plot(numRange,power2M,'ro-',numRange,power1M,'bx-',numRange,power2L,'r.:',numRange,power1L,'b.:',numRange,power3,'g.-',numRange,power4,'c.-','LineWidth',1);
    end
    xlim([numRange(1) numRange(end)]);
    ylim([0 1]);
    title(titlechar);
end
suptitle('Independence Test Powers for 20 Simulated Dependencies with Increasing Sample Size')
% xlabel('Sample Size');
% ylabel('Empirical Testing Power');
switch option
    case 0
        legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Location','SouthOutside');
    case 1
        legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','SouthOutside');
end

figure
s=5;
t=4;
for j=1:total
    filename=strcat('CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=' Data';
    switch type
        case 1
            titlechar=' 1. Linear';
        case 2
            titlechar=' 2. Quadratic';
        case 3
            titlechar=' 3. Cubic';
        case 4
            titlechar=' 4. Sine Period 1/2';
        case 5
            titlechar=' 5. Sine Period 1/8';
        case 6
            titlechar=' 6. X\^(1/4)';
        case 7
            titlechar=' 7. Circle';
        case 8
            titlechar=' 8. Step Function';
        case 9
            titlechar=' 9. Exp(X)';
        case 10
            titlechar=' 10. Uncorrelated Binomial';
        case 11
            titlechar=' 11. W';
        case 12
            titlechar=' 12. Square';
        case 13
            titlechar=' 13. Diamond';
        case 14
            titlechar=' 14. Parabola';
        case 15
            titlechar=' 15. Two Parabolas';
        case 16
            titlechar=' 16. Circle 2';
        case 17
            titlechar=' 17. Independent Clouds';
        case 18
            titlechar=' 18. Joint Normal';
        case 19
            titlechar=' 19. Log(X^2)';
        case 20
            titlechar=' 20. Multiplicative Noise';
    end
    %     %%Plot the power w.r.t. neighborhood
    for i=1:length(numRange)
        p1(i,:)=max(power1(:,:,i),[],1);
        p2(i,:)=max(power2(:,:,i),[],1);
        %         p1(i,:)=power1(i,numRange(i),:);%
        %         p2(i,:)=power2(i,numRange(i),:);%
        %         p3(i,:)=power3(i,numRange(i),:);%check
        power1L(i)=power1(numRange(i),numRange(i),i);
        power2L(i)=power2(numRange(i),numRange(i),i);
    end
    power1=p1;power2=p2;
    K=n;kmin=1;thres=0.8;
    ind=[find(max(power1,[],2)>thres,1) find(max(power2,[],2)>thres,1) find(power3>thres,1) find(power4>thres,1) lim];
    lim=min(ind);
    xaxis=kmin:numRange(lim);
    if option==0
        plot(xaxis,power2(lim,xaxis),'r.-',xaxis, power1(lim,xaxis),'b.-',xaxis,power2L(lim)*ones(length(xaxis),1),'r.:',xaxis, power1L(lim)*ones(length(xaxis),1),'b.:',xaxis,power3(lim)*ones(length(xaxis),1),'g.-','LineWidth',1);
    else
        plot(xaxis,power2(lim,xaxis),'r.-',xaxis, power1(lim,xaxis),'b.-',xaxis,power2L(lim)*ones(length(xaxis),1),'r.:',xaxis, power1L(lim)*ones(length(xaxis),1),'b.:',xaxis,power3(lim)*ones(length(xaxis),1),'g.-',xaxis,power4(lim)*ones(length(xaxis),1),'c.-','LineWidth',1);
    end
    xlim([kmin numRange(lim)]);
    ylim([0 1]);
    title(titlechar);
end
suptitle('Independence Test Powers for 20 Simulated Dependencies with Increasing Neighborhood at Fixed Sample Size')
% xlabel('Neighborhood Size');
% ylabel('Empirical Testing Power');
if option==0
    legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Location','SouthOutside');
else
    legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','SouthOutside');   
end

%%%performance profile
figure
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
profile=zeros(6,length(xaxis));
%load data
for j=1:total
    filename=strcat('CorrIndTestType',num2str(j),'N100Dim1.mat');
    load(filename)
    for i=1:length(numRange)
        power1L(i)=power1(numRange(i),numRange(i),i);
        power2L(i)=power2(numRange(i),numRange(i),i);
        power1M(i)=max(max(power1(:,:,i)));
        power2M(i)=max(max(power2(:,:,i)));
    end
    thres=0.8;
    ind=[find(power1M>thres,1) find(power2M>thres,1) find(power3>thres,1) find(power4>thres,1) lim];
    lim=min(ind);
    pos=lim;
    switch option
        case 0
            power=[power2M(pos), power1M(pos), power2L(pos), power1L(pos),power3(pos)];
        case 1
            power=[power2M(pos), power1M(pos), power2L(pos), power1L(pos),power3(pos),power4(pos)];
    end
    pmax=max(power);
    for k=1:length(power)
        tmp=(pmax-power(k));
        tmpInd=ceil(tmp/interval)+1;
        profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
    end
end
%
profile=profile./total;
sumP=ceil(mean(profile,2)*1000)/1000;
switch option
    case 0
        plot(xaxis,profile(1,:),'r.-',xaxis, profile(2,:),'b.-',xaxis,profile(3,:),'r.:',xaxis, profile(4,:),'b.:',xaxis,profile(5,:),'g.-','LineWidth',2);
        legend(strcat('Local Modified Distance Correlation, AUC=', num2str(sumP(1))),strcat('Local Original Distance Correlation, AUC=', num2str(sumP(2))),strcat('Modified Distance Correlation, AUC=', num2str(sumP(3))),strcat('Original Distance Correlation, AUC=', num2str(sumP(4))),strcat('HHG, AUC=', num2str(sumP(5))),'Location','SouthEast');
    case 1
        plot(xaxis,profile(1,:),'r.-',xaxis, profile(2,:),'b.-',xaxis,profile(3,:),'r.:',xaxis, profile(4,:),'b.:',xaxis,profile(5,:),'g.-',xaxis,profile(6,:),'c.-','LineWidth',2);
        legend(strcat('Local Modified Distance Correlation, AUC=', num2str(sumP(1))),strcat('Local Original Distance Correlation, AUC=', num2str(sumP(2))),strcat('Modified Distance Correlation, AUC=', num2str(sumP(3))),strcat('Original Distance Correlation, AUC=', num2str(sumP(4))),strcat('HHG, AUC=', num2str(sumP(5))),strcat('Mantel, AUC=', num2str(sumP(6))),'Location','SouthEast');
end
xlabel('Difference with the Best Method','FontSize',12);
ylabel('Performance Profiles','FontSize',12);
ylim([0 1]);
titleStr = strcat('Performance Profiles for Dimension 1');
title(titleStr,'FontSize',12);
%

%%%performance profile AUC for n
figure
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
limN=10;
sumP=zeros(6,limN);
for ll=1:limN
    profile=zeros(6,length(xaxis));
    %load data
    for j=1:total
        filename=strcat('CorrIndTestType',num2str(j),'N100Dim1.mat');
        load(filename)
        for i=1:length(numRange)
            power1L(i)=power1(numRange(i),numRange(i),i);
            power2L(i)=power2(numRange(i),numRange(i),i);
            power1M(i)=max(max(power1(:,:,i)));
            power2M(i)=max(max(power2(:,:,i)));
        end
        thres=ll/limN;%thres=0.8;
        ind=[find(power1M>thres,1) find(power2M>thres,1) find(power3>thres,1) find(power4>thres,1) lim];
        lim=min(ind);
        pos=lim;
        switch option
            case 0
                power=[power2M(pos), power1M(pos), power2L(pos), power1L(pos),power3(pos)];
            case 1
                power=[power2M(pos), power1M(pos), power2L(pos), power1L(pos),power3(pos),power4(pos)];
        end
        pmax=max(power);
        for k=1:length(power)
            tmp=(pmax-power(k));
            tmpInd=ceil(tmp/interval)+1;
            profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
        end
    end
    %
    profile=profile./total;
    sumP(:,ll)=ceil(mean(profile,2)*1000)/1000;
end
xaxis=1/limN:1/limN:1;
%xaxis=numRange;
switch option
    case 0
        plot(xaxis,sumP(1,:),'r.-',xaxis, sumP(2,:),'b.-',xaxis,sumP(3,:),'r.:',xaxis, sumP(4,:),'b.:',xaxis,sumP(5,:),'g.-','LineWidth',2);
        legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Location','SouthEast');
    case 1
        plot(xaxis,sumP(1,:),'r.-',xaxis, sumP(2,:),'b.-',xaxis,sumP(3,:),'r.:',xaxis, sumP(4,:),'b.:',xaxis,sumP(5,:),'g.-',xaxis,sumP(6,:),'c.-','LineWidth',2);
        legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','SouthEast');
end
%add rdcorr at k=n-1
%plot(xaxis,sumP(1,:),'r.-',xaxis, sumP(2,:),'b.-',xaxis,sumP(3,:),'g.-',xaxis,sumP(4,:),'g.-',xaxis,sumP(5,:),'m.-','LineWidth',2);
%legend('Rank Distance Correlation','Original Distance Correlation','Modified Distance Correlation','HHG','Location','SouthEast');% Figure title/labels
xlabel('Threshold','FontSize',12);
%xlabel('Sample Size');
ylabel('AUC','FontSize',12);
%xlim([numRange(1) numRange(end)]);
ylim([0 1]);
titleStr = strcat('Performance Profiles AUC for Dimension 1');
title(titleStr,'FontSize',12);
%


%Plot 5-8
figure
s=5;
t=4;
for j=1:total
    filename=strcat('CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=' Data';
    switch type
        case 1
            titlechar=' 1. Linear';
        case 2
            titlechar=' 2. Quadratic';
        case 3
            titlechar=' 3. Cubic';
        case 4
            titlechar=' 4. Sine Period 1/2';
        case 5
            titlechar=' 5. Sine Period 1/8';
        case 6
            titlechar=' 6. X\^(1/4)';
        case 7
            titlechar=' 7. Circle';
        case 8
            titlechar=' 8. Step Function';
        case 9
            titlechar=' 9. Exp(X)';
        case 10
            titlechar=' 10. Uncorrelated Binomial';
        case 11
            titlechar=' 11. W';
        case 12
            titlechar=' 12. Square';
        case 13
            titlechar=' 13. Diamond';
        case 14
            titlechar=' 14. Parabola';
        case 15
            titlechar=' 15. Two Parabolas';
        case 16
            titlechar=' 16. Circle 2';
        case 17
            titlechar=' 17. Independent Clouds';
        case 18
            titlechar=' 18. Joint Normal';
        case 19
            titlechar=' 19. Log(X^2)';
        case 20
            titlechar=' 20. Multiplicative Noise';
    end
    power1M=zeros(length(dimRange),1);power2M=zeros(length(dimRange),1);power1L=zeros(length(dimRange),1);power2L=zeros(length(dimRange),1);
    for i=1:length(dimRange)
        power1L(i)=power1(n,n,i);
        power2L(i)=power2(n,n,i);
        power1M(i)=max(max(power1(:,:,i)));
        power2M(i)=max(max(power2(:,:,i)));
    end
    switch option
        case 0
            plot(dimRange,power2M,'ro-',dimRange,power1M,'bx-',dimRange,power2L,'r.:',dimRange,power1L,'b.:',dimRange,power3,'g.-','LineWidth',1);
        case 1
            plot(dimRange,power2M,'ro-',dimRange,power1M,'bx-',dimRange,power2L,'r.:',dimRange,power1L,'b.:',dimRange,power3,'g.-',dimRange,power4,'c.-','LineWidth',1);
    end
    xlim([dimRange(1) dimRange(end)]);
    ylim([0 1]);
    title(titlechar);
end
suptitle('Independence Test Powers for 20 Simulated Dependencies with Increasing Dimension')
% xlabel('Sample Size');
% ylabel('Empirical Testing Power');
switch option
    case 0
        legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Location','SouthOutside');
    case 1
        legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','SouthOutside');
end

figure
s=5;
t=4;
for j=1:total
    filename=strcat('CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    subplot(s,t,j)
    titlechar=' Data';
    switch type
        case 1
            titlechar=' 1. Linear';
        case 2
            titlechar=' 2. Quadratic';
        case 3
            titlechar=' 3. Cubic';
        case 4
            titlechar=' 4. Sine Period 1/2';
        case 5
            titlechar=' 5. Sine Period 1/8';
        case 6
            titlechar=' 6. X\^(1/4)';
        case 7
            titlechar=' 7. Circle';
        case 8
            titlechar=' 8. Step Function';
        case 9
            titlechar=' 9. Exp(X)';
        case 10
            titlechar=' 10. Uncorrelated Binomial';
        case 11
            titlechar=' 11. W';
        case 12
            titlechar=' 12. Square';
        case 13
            titlechar=' 13. Diamond';
        case 14
            titlechar=' 14. Parabola';
        case 15
            titlechar=' 15. Two Parabolas';
        case 16
            titlechar=' 16. Circle 2';
        case 17
            titlechar=' 17. Independent Clouds';
        case 18
            titlechar=' 18. Joint Normal';
        case 19
            titlechar=' 19. Log(X^2)';
        case 20
            titlechar=' 20. Multiplicative Noise';
    end
    %     %%Plot the power w.r.t. neighborhood
    for i=1:length(dimRange)
        p1(i,:)=max(power1(:,:,i),[],1);
        p2(i,:)=max(power2(:,:,i),[],1);
        power1L(i)=power1(end,end,i);
        power2L(i)=power2(end,end,i);
    end
    power1=p1;power2=p2;
    K=n;kmin=1;thres=0.5;
    ind=[find(max(power2,[],2)==max(max(power2)),1,'first') find(max(power2,[],2)>thres,1,'last') find(max(power1,[],2)>thres,1,'last') find(max(power3,[],2)>thres,1,'last') find(max(power4,[],2)>thres,1,'last')];
    lim=max(ind);
    xaxis=kmin:K;
    if option==0;
        plot(xaxis,power2(lim,xaxis),'r.-',xaxis, power1(lim,xaxis),'b.-',xaxis,power2L(lim)*ones(length(xaxis),1),'r.:',xaxis,power1L(lim)*ones(length(xaxis),1),'b.:',xaxis,power3(lim)*ones(length(xaxis),1),'g.-','LineWidth',1);
    else
        plot(xaxis,power2(lim,xaxis),'r.-',xaxis, power1(lim,xaxis),'b.-',xaxis,power2L(lim)*ones(length(xaxis),1),'r.:',xaxis,power1L(lim)*ones(length(xaxis),1),'b.:',xaxis,power3(lim)*ones(length(xaxis),1),'g.-',xaxis,power4(lim)*ones(length(xaxis),1),'c.-','LineWidth',1);
    end
    xlim([kmin K]);
    ylim([0 1]);
    title(titlechar);
end
suptitle('Independence Test Powers for 20 Simulated Dependencies with Increasing Neighborhood at Fixed Dimension')
% xlabel('Neighborhood Size');
% ylabel('Empirical Testing Power');
if option==0
    legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Location','SouthOutside');
else
    legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','SouthOutside');
end

%%%performance profile
figure
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
profile=zeros(6,length(xaxis));
%load data
for j=1:total
    filename=strcat('CorrIndTestDimType',num2str(j),'N100Dim.mat');
    load(filename)
    thres=0.5;
    for i=1:length(dimRange)
        power1L(i)=power1(n,n,i);
        power2L(i)=power2(n,n,i);
        power1M(i)=max(max(power1(:,:,i)));
        power2M(i)=max(max(power2(:,:,i)));
    end
    ind=[find(power1M>thres,1,'last') find(power2M>thres,1,'last') find(power3>thres,1,'last') find(power4>thres,1,'last') 1];
    lim=max(ind);
    pos=lim;
    switch option
        case 0
            power=[power2M(pos), power1M(pos), power2L(pos), power1L(pos),power3(pos)];
        case 1
            power=[power2M(pos), power1M(pos), power2L(pos), power1L(pos),power3(pos),power4(pos)];
    end
    pmax=max(power);
    for k=1:length(power)
        tmp=(pmax-power(k));
        tmpInd=ceil(tmp/interval)+1;
        profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
    end
end
%
profile=profile./total;
sumP=ceil(mean(profile,2)*1000)/1000;
switch option
    case 0
        plot(xaxis,profile(1,:),'r.-',xaxis, profile(2,:),'b.-',xaxis,profile(3,:),'r.:',xaxis, profile(4,:),'b.:',xaxis,profile(5,:),'g.-','LineWidth',2);
        legend(strcat('Local Modified Distance Correlation, AUC=', num2str(sumP(1))),strcat('Local Original Distance Correlation, AUC=', num2str(sumP(2))),strcat('Modified Distance Correlation, AUC=', num2str(sumP(3))),strcat('Original Distance Correlation, AUC=', num2str(sumP(4))),strcat('HHG, AUC=', num2str(sumP(5))),'Location','SouthEast');
    case 1
        plot(xaxis,profile(1,:),'r.-',xaxis, profile(2,:),'b.-',xaxis,profile(3,:),'r.:',xaxis, profile(4,:),'b.:',xaxis,profile(5,:),'g.-',xaxis,profile(6,:),'c.-','LineWidth',2);
        legend(strcat('Local Modified Distance Correlation, AUC=', num2str(sumP(1))),strcat('Local Original Distance Correlation, AUC=', num2str(sumP(2))),strcat('Modified Distance Correlation, AUC=', num2str(sumP(3))),strcat('Original Distance Correlation, AUC=', num2str(sumP(4))),strcat('HHG, AUC=', num2str(sumP(5))),strcat('Mantel, AUC=', num2str(sumP(6))),'Location','SouthEast');
end
%add rdcorr at k=n-1
%plot(xaxis,profile(1,:),'r.-',xaxis, profile(2,:),'b.-',xaxis,profile(3,:),'g.-',xaxis,profile(4,:),'g.-',xaxis,profile(5,:),'m.-','LineWidth',2);
%legend(strcat('Rank Distance Correlation*, AUC=', num2str(sumP(1))),strcat('Original Distance Correlation*, AUC=', num2str(sumP(2))),strcat('Modified Distance Correlation*, AUC=', num2str(sumP(3))),strcat('HHG, AUC=', num2str(sumP(4))),strcat('Rank Distance Correlation at n-1, AUC=', num2str(sumP(5))), 'Location','SouthEast');
% Figure title/labels
xlabel('Difference with the Best Method','FontSize',12);
ylabel('Performance Profiles','FontSize',12);
ylim([0 1]);
titleStr = strcat('Performance Profiles for Increasing Dimension');
title(titleStr,'FontSize',12);

%%%performance profile
figure
a=0;b=1;interval=0.01;
xaxis=a:interval:b;
limN=10;
sumP=zeros(6,limN);
for ll=1:limN
    profile=zeros(6,length(xaxis));
    %load data
    for j=1:total
        filename=strcat('CorrIndTestDimType',num2str(j),'N100Dim.mat');
        load(filename)
        thres=ll/limN;
        for i=1:length(dimRange)
            power1L(i)=power1(n,n,i);
            power2L(i)=power2(n,n,i);
            power1M(i)=max(max(power1(:,:,i)));
            power2M(i)=max(max(power2(:,:,i)));
        end
        ind=[find(power1M>thres,1,'last') find(power2M>thres,1,'last') find(power3>thres,1,'last') find(power4>thres,1,'last') 1];
        lim=max(ind);
        pos=lim;
        switch option
            case 0
                power=[power2M(pos), power1M(pos), power2L(pos), power1L(pos),power3(pos)];
            case 1
                power=[power2M(pos), power1M(pos), power2L(pos), power1L(pos),power3(pos),power4(pos)];
        end
        pmax=max(power);
        for k=1:length(power)
            tmp=(pmax-power(k));
            tmpInd=ceil(tmp/interval)+1;
            profile(k,tmpInd:end)=profile(k,tmpInd:end)+1;
        end
    end
    %
    profile=profile./total;
    sumP(:,ll)=ceil(mean(profile,2)*1000)/1000;
end
xaxis=1/limN:1/limN:1;
switch option
    case 0
        plot(xaxis,sumP(1,:),'r.-',xaxis, sumP(2,:),'b.-',xaxis,sumP(3,:),'r.:',xaxis, sumP(4,:),'b.:',xaxis,sumP(5,:),'g.-','LineWidth',2);
        legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Location','SouthEast');
    case 1
        plot(xaxis,sumP(1,:),'r.-',xaxis, sumP(2,:),'b.-',xaxis,sumP(3,:),'r.:',xaxis, sumP(4,:),'b.:',xaxis,sumP(5,:),'g.-',xaxis,sumP(6,:),'c.-','LineWidth',2);
        legend('Local Modified Distance Correlation','Local Original Distance Correlation','Modified Distance Correlation','Original Distance Correlation','HHG','Mantel','Location','SouthEast');
end
%add rdcorr at k=n-1
%plot(xaxis,sumP(1,:),'r.-',xaxis, sumP(2,:),'b.-',xaxis,sumP(3,:),'g.-',xaxis,sumP(4,:),'g.-',xaxis,sumP(5,:),'m.-','LineWidth',2);
%legend('Rank Distance Correlation*','Original Distance Correlation*','Modified Distance Correlation*','HHG','Location','SouthEast');% Figure title/labels
xlabel('Threshold','FontSize',12);
ylabel('AUC','FontSize',12);
ylim([0 1]);
titleStr = strcat('Performance Profiles AUC for Increasing Dimension');
title(titleStr,'FontSize',12);
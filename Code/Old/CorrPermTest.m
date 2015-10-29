function [p1, p2, p3, p4]=CorrPermTest(type,n,dim,lim,rep1,rep2,noise)
% Author: Cencheng Shen
% Permutation Tests for identifying correlation, returning p-value with respect to increasing sample size.
% The output is the p-value of rankDCorr, dCorr, modified dCorr, HHG.
%
% Type 1-10: for Y=f(XA), similar to Tibs comment 2008
% n=100; dim=100; lim=10; rep1=100; rep2=1000;
% CorrPermTest(1,n,dim,lim,rep1,rep2);
%
% Type 11-17: seven correlation types in Wikipedia correlation figure /
% Table 3 of HHG, which uses dim=1.
% n=100; dim=1; lim=10; rep1=100; rep2=1000;
% CorrPermTest(11,n,dim,lim,rep1,rep2);
%
% Type 18-20: 3 correlation types in Szeley 2007, which uses dim=5.
% n=100; dim=5; lim=10; rep1=100; rep2=1000;
% CorrPermTest(18,n,dim,lim,rep1,rep2);
%
% Input real data into type for: BrainCxP
% load('BrainCP')
% dim=5;
% x=embeddingMethod(distC, 'Dis', dim);
% y=embeddingMethod(distP, 'Dis', dim);
% n=42; dim=5; lim=7; rep1=100; rep2=1000;
% CorrPermTest([x y],n,dim,lim,rep1,rep2);

% Parameters:
% In the input, rep1 specifies the number of MC replicates.
% rep2 specifies the number of random permutations.
display=0;
if nargin<7
    noise=1; %noise level default
end
K=n-1;
alpha=0.05; %type 1 error level
option1=1; option2=1; option3=1; option4=1; %Control whether to calculate the respective correlation statistic or not.

if lim==0
    numRange=n; %do sample size at n
else
    numRange=ceil(n/lim):ceil(n/lim):n; %do sample size at the interval of ceil(n/lim).
end
lim=length(numRange);

% High-dimensional decay
A=ones(dim,1);
%A=A./(ceil(dim*rand(dim,1))); %random decay
for d=1:dim
    A(d)=A(d)/d; %fixed decay
end

% Output
p1=zeros(lim,K,rep1); p2=zeros(lim,K,rep1);p3=zeros(lim,K,rep1);p4=zeros(lim,rep1);%p-values for rankdCorr, dCorr, modified dCorr, HHG
power1=zeros(lim,K,1); power2=zeros(lim,K,1);power3=zeros(lim,K,1);power4=zeros(lim,1);%average powers for rankdCorr, dCorr, modified dCorr, HHG
dCor1=zeros(rep2,K);dCor2=zeros(rep2,K);dCor3=zeros(rep2,K);dCor4=zeros(rep2,1);
d=dim;

%rep1 MC replicates
for r1=1:rep1
    % Add noise if necessary; by default noiseless when noise=0.
    if noise ~=0
        eps=mvnrnd(0,1,n); %Gaussian noise
    else
        eps=zeros(n,1);
    end
    
    if size(type,1)>1 %Real data testing
        x=type(:,1:d);
        y=type(:,d+1:2*d);
    else
        x=unifrnd(-1,1,n,d);
        xA=x*A;
        switch type %Various simulations and their example parameters
            case 1 %Linear, dim=100, n=100
                y=xA+1*noise*eps;
            case 2 %Quadratic, dim=100, n=200
                y=4*(xA-0.5).^2+4*noise*eps;
            case 3 %Cubic, dim=100, n=200
                y=128*(xA-1/3).^3+48*(xA-1/3).^2-12*(xA-1/3)+80*noise*eps;
            case 4 %Sine 1/2, dim=1, n=100
                y=sin(4*pi*xA/max(xA))+2*noise*eps;
            case 5 %Sine 1/8, dim=1, n=200
                y=sin(16*pi*xA/max(xA))+noise*eps;
                %y=(xA-1/3).^10+1000*noise*eps;
            case 6 %X^0.25, dim=1, n=40
                y=abs(xA).^0.25+noise/5*eps;
            case 7 %Circle, dim=10, n=100
                y=(binornd(1,0.5,n,1)*2-1).*abs(1-(xA/max(xA)*2-1).^2).^0.5+noise*eps;
            case 8 %Step Function, dim=100, n=100
                y=(xA>mean(xA))+1*noise*eps;
            case 9 %Exponential, dim=100, n=100
                y=exp(xA)+1.5*noise*eps;
            case 10 %Uncorrelated X & Y, dim=1, n=40
                x=binornd(1,0.5,n,d);
                y=(binornd(1,0.5,n,d)*2-1)*A;
                y=x*A.*y+0.6*noise*eps;
            case 11 %W Shape, dim=1, n=100
                x=x+unifrnd(0,1,n,d)/3;
                y=4*( ( xA.^2 - 1/2 ).^2 + unifrnd(0,1,n,d)*A/500 )+0.3*noise*eps;
            case 12 %Square, dim=1, n=200
                u=unifrnd(-1,1,n,d)*A;
                v=unifrnd(-1,1,n,d)*A;
                theta=-pi/8;
                tmp=[cos(theta) -sin(theta); sin(theta) cos(theta)];
                uv=[u v] * tmp;
                x=uv(:,1);
                y=uv(:,2)+noise*eps;
            case 13 %Diamond, dim=1, n=200
                u=unifrnd(-1,1,n,d)*A;
                v=unifrnd(-1,1,n,d)*A;
                theta=-pi/4;
                tmp=[cos(theta) -sin(theta); sin(theta) cos(theta)];
                uv=[u v] * tmp;
                x=uv(:,1);
                y=uv(:,2)+noise*eps;
            case 14 %Parabola, dim=1, n=100
                y=( xA.^2  + unifrnd(0,1,n,1))/2+0.2*noise*eps;
            case 15 %Two Parabolas, dim=1, n=100
                y=( xA.^2  + unifrnd(0,1,n,1))/2.*(randsample(2, n, true)-1.5)*2+noise*eps;
            case 16 %Circle, dim=1, n=100
                x=sin(xA*pi)+mvnrnd(0,1,n)/8;
                y=cos(xA*pi)+mvnrnd(0,1,n)/8+0.2*noise*eps;
            case 17 %Independent clouds, dim=1, n=100
                x=mvnrnd(zeros(n,d),eye(d),n)*A/3+(randsample(2, n, true)-1.5)*2;
                y=mvnrnd(zeros(n,d),eye(d),n)*A/3+(randsample(2, n, true)-1.5)*2+noise*eps;
            case 18 %Joint Normal, dim=5, n=200
                rho=1/(dim*2);
                cov1=[eye(d) rho*ones(d)];
                cov2=[rho*ones(d) eye(d)];
                %                     cov1=[eye(d) rho*eye(d)];
                %                     cov2=[rho*eye(d) eye(d)];
                cov=[cov1' cov2'];
                x=mvnrnd(zeros(n,2*d),cov,n);
                y=x(:,d+1:2*d)+2*noise*repmat(eps,1,d);
                x=x(:,1:d);
            case 19 %Log(X^2), dim=5, n=100
                x=mvnrnd(zeros(n, d),eye(d));
                y=log(x.^2)+2*noise*repmat(eps,1,d);
            case 20 %Multiplicative Noise, dim=5, n=100
                x=mvnrnd(zeros(n, d),eye(d));
                y=mvnrnd(zeros(n, 1),eye(1));
                y=x.*repmat(y,1,d)+2*noise*repmat(eps,1,d);
        end
    end
    
    % Iterate through all sample size choices
    C1=squareform(pdist(x));
    P1=squareform(pdist(y));
    for i=1:lim
        nn=numRange(i);
        C=C1(1:nn,1:nn);
        P=P1(1:nn,1:nn);
        disRankC=disToRanks(C);
        disRankP=disToRanks(P);
        
        % Permute the second dataset for rep2 times, and calculate the 4 correlation statistics
        for r2=1:rep2
            per=randperm(nn);
            Pa=P(per,per);
            disRank=[disRankC disRankP(per, per)];
            if option1~=0
                dCor1(r2,1:nn-1)=rankDCorr(C,Pa,disRank,1);
            end
            if option2~=0
                dCor2(r2,1:nn-1)=rankDCorr(C,Pa,disRank,2);
            end
            if option3~=0
                dCor3(r2,1:nn-1)=rankDCorr(C,Pa,disRank,3);
            end
            if option4~=0
                dCor4(r2)=HHG(C,Pa);
            end
        end
        
        % Calculate the correlation statistics for the original data, and
        % derive the p-value for the permutation test.
        disRank=[disRankC disRankP];
        if option1~=0
            cut1=rankDCorr(C,P,disRank,1);
            p1(i,1:nn-1,r1)=mean(dCor1(:,1:nn-1)<repmat(cut1,1,rep2)',1);
        end
        if option2~=0
            cut2=rankDCorr(C,P,disRank,2);
            p2(i,1:nn-1,r1)=mean(dCor2(:,1:nn-1)<repmat(cut2,1,rep2)',1);
        end
        if option3~=0
            cut3=rankDCorr(C,P,disRank,3);
            p3(i,1:nn-1,r1)=mean(dCor3(:,1:nn-1)<repmat(cut3,1,rep2)',1);
        end
        if option4~=0
            cut4=HHG(C,P);
            p4(i,r1)=mean(dCor4<cut4);
        end
    end
    power1=power1+(p1(:,:,r1)>(1-alpha))/rep1;
    power2=power2+(p2(:,:,r1)>(1-alpha))/rep1;
    power3=power3+(p3(:,:,r1)>(1-alpha))/rep1;
    power4=power4+(p4(:,r1)>(1-alpha))/rep1;
end

% Output the mean p-value and standard deviation. Std is not meaningful
% when rep1=1.
std1=std(p1,0,3);
std2=std(p2,0,3);
std3=std(p3,0,3);
std4=std(p4,0,2);
p1=1-mean(p1,3);
p2=1-mean(p2,3);
p3=1-mean(p3,3);
p4=1-mean(p4,2);

% Save the results
if size(type,1)>1
    type=' Brain CxP';
end
filename=strcat('CorrPermTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
save(filename,'power1','power2','power3','power4','p1','p2','p3','p4','std1','std2','std3','std4','dCor1','dCor2','dCor3','dCor4','cut1','cut2','cut3','cut4','type','n','numRange','dim','noise','lim','rep1','rep2');
%save(filename,'power1','power2','power3','power4','p1','p2','p3','p4','std1','std2','std3','std4','type','n','numRange','dim','noise','lim','rep1','rep2');

% Display and save picture. By default do not display.
if display~=0
    filename=strcat('CorrPermTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
    titlechar=type;
    switch type
        case 1
            titlechar=' Linear';
        case 2
            titlechar=' Quadratic';
        case 3
            titlechar=' Cubic';
        case 4
            titlechar=' Sine Period 1/2';
        case 5
            titlechar=' Sine Period 1/8';
        case 6
            titlechar=' X\^(1/4)';
        case 7
            titlechar=' Circle';
        case 8
            titlechar=' Step Function';
        case 9
            titlechar=' Exp(X)';
        case 10
            titlechar=' Uncorrelated Binomial';
        case 11
            titlechar=' W';
        case 12
            titlechar=' Square';
        case 13
            titlechar=' Diamond';
        case 14
            titlechar=' Parabola';
        case 15
            titlechar=' Two Parabolas';
        case 16
            titlechar=' Circle';
        case 17
            titlechar=' Independent Clouds';
        case 18
            titlechar=' Joint Normal';
        case 19
            titlechar=' Log(X^2)';
        case 20
            titlechar=' Multiplicative Noise';
    end

    figure
    %%% Plot the average p-value
    % plot(numRange,min(p1(:,:),[],2),'ro-',numRange, min(p2(:,:),[],2),'bx-',numRange,min(p3(:,:),[],2),'g^-',numRange,p4,'c.-','LineWidth',2);
    % plot(numRange,min(p1(:,:),[],2),'ro-',numRange, p2(:,end),'bx-',numRange, p3(:,end),'g^-',numRange,p4,'c.-','LineWidth',2);
    %legend('Rank Distance Correlation','Distance Correlation','Modified Distance Correlation','HHG','Location','NorthEast');
    % Plot the permutation test power
    %plot(numRange,max(power1,[],2),'ro-',numRange, max(power2,[],2),'bx-',numRange,max(power3,[],2),'g^-',numRange,power4,'c.-','LineWidth',2);
    for i=1:length(numRange)
        power2L(i)=power2(i,numRange(i)-1);
        power3L(i)=power3(i,numRange(i)-1);
    end
    plot(numRange,max(power1,[],2),'ro-',numRange, power2L,'bx-',numRange, power3L,'g^-',numRange,power4,'c.-','LineWidth',2);
    legend('Rank Distance Correlation','Distance Correlation','Modified Distance Correlation','HHG','Location','SouthEast');
    xlim([numRange(1) numRange(end)]);
    ylim([0 1]);
    xlabel('Sample Size');
    ylabel('Power');
    titleStr = strcat('Permutation Test for ', titlechar,' upto n=', num2str(n));
    title(titleStr);
    saveas(gcf,filename,'jpeg');
    
    %     %%Plot the power/p-value w.r.t. neighborhood
    %     figure
    %     K=n-1;
    %     kmin=ceil(K/10);xaxis=kmin:K;
    %     %%plot(xaxis,power1(lim,xaxis),'ro-',xaxis, power2(lim)*ones(length(xaxis),1),'bx-',xaxis,power3(lim)*ones(length(xaxis),1),'g^-',xaxis,power4(lim)*ones(length(xaxis),1),'c.-','LineWidth',1);
    %     plot(xaxis,p1(lim,xaxis),'ro-',xaxis, p2(lim)*ones(length(xaxis),1),'bx-',xaxis,p3(lim)*ones(length(xaxis),1),'g^-',xaxis,p4(lim)*ones(length(xaxis),1),'c.-','LineWidth',1);
    %     xlabel('Neighborhood Size');
    %     xlim([1 K]);
    %     ylabel('P-Value');
    %     ylim([0 1]);
    
    %     % Figure title/labels
    %     legend('Rank Distance Correlation','Distance Correlation','Modified Distance Correlation','HHG','Location','NorthEast');
    %     titleStr = strcat('Permutation Test for ', titlechar,' at n=', num2str(n));
    %     title(titleStr);
    %     saveas(gcf,filename,'jpeg');
end
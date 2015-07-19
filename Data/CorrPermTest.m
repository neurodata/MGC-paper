function [power1, power2, power3, power4]=CorrPermTest(type,n,dim,lim,rep1,rep2)
% Author: Cencheng Shen
% Permutation Tests for identifying correlation, returning p-value with respect to increasing sample size.
% The output is the p-value of dCorr, rankDCorr, modified dCorr, HHG.
%
% Type 1-10: for Y=f(XA), similar to Tibs comment 2008
% n=100; dim=100; lim=10; rep1=1; rep2=1000;
% CorrPermTest(1,n,dim,lim,rep1,rep2);
%
% Type 11-17: seven correlation types in Wikipedia correlation figure /
% Table 3 of HHG, which uses dim=1.
% n=100; dim=1; lim=10; rep1=1; rep2=1000;
% CorrPermTest(11,n,dim,lim,rep1,rep2);
%
% Type 18-20: 3 correlation types in Szeley 2007, which uses dim=5.
% n=100; dim=5; lim=10; rep1=1; rep2=1000;
% CorrPermTest(18,n,dim,lim,rep1,rep2);
%
% Input real data into type for: BrainCxP
% load('BrainCP')
% dim=5;
% x=embeddingMethod(distC, 'Dis', dim);
% y=embeddingMethod(distP, 'Dis', dim);
% n=42; dim=5; lim=7; rep1=100; rep2=1000;
% CorrPermTest([x y],n,dim,lim,rep1,rep2);
%
% CorrPermTest([X_data(:,1:5000), [Y_data(:,1:5000)', zeros(5000,1)]'],100, 3,5);
% CorrPermTest([squareform(pdist(x)) squareform(pdist(y))],n, dim, n-1);

% Parameters:
% In the input, rep1 specifies the number of MC replicates.
% rep2 specifies the number of random permutations.
display=0;
noise=0;
K=n-1;

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
power1=zeros(lim,rep1); power2=zeros(lim,K,rep1);power3=zeros(lim,rep1);power4=zeros(lim,rep1);%p-values for dCorr, rankdCorr, modified dCorr, HHG
dCor1=zeros(rep2,1);dCor2=zeros(rep2,K);dCor3=zeros(rep2,1);dCor4=zeros(rep2,1);
d=dim;

for r1=1:rep1
    % Add noise if necessary; by default noiseless when noise=0.
    if noise ~=0
        eps=mvnrnd(0,1,n); %Gaussian noise
    else
        eps=0;
    end
    
    if size(type,1)>1 %Real data testing
        x=type(:,1:d);
        y=type(:,d+1:2*d);
    else
        x=unifrnd(-1,1,n,d);
        xA=x*A;
        switch type %Various simulations and their example parameters
            case 1 %Linear, dim=100, n=100
                y=xA+noise*eps;
            case 2 %Quadratic, dim=100, n=200
                y=4*(xA-0.5).^2+10*noise*eps;
            case 3 %Cubic, dim=100, n=200
                y=128*(xA-1/3).^3+48*(xA-1/3).^2-12*(xA-1/3)+1000*noise*eps;
            case 4 %Sine 1/2, dim=1, n=100
                y=sin(4*pi*xA/max(xA))+2*noise*eps;
            case 5 %Sine 1/8, dim=1, n=200
                y=sin(16*pi*xA/max(xA))+noise*eps;
                %y=(xA-1/3).^10+1000*noise*eps;
            case 6 %X^0.25, dim=1, n=50
                y=abs(xA).^0.25+noise*eps;
            case 7 %Circle, dim=10, n=100
                y=(binornd(1,0.5,n,1)*2-1).*abs(1-(xA/max(xA)*2-1).^2).^0.5+noise/4*eps;
            case 8 %Step Function, dim=100, n=100
                y=(xA>mean(xA))+1*noise*eps;
            case 9 %Exponential, dim=100, n=100
                y=exp(xA)+10*noise*eps;
            case 10 %Log(X), dim=1, n=50
                y=log(abs(xA)+1)+noise*eps;
            case 11 %W Shape, dim=1, n=100
                x=x+unifrnd(0,1,n,d)/3;
                y=4*( ( xA.^2 - 1/2 ).^2 + unifrnd(0,1,n,d)*A/500 )+noise*eps;
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
                y=( xA.^2  + unifrnd(0,1,n,1))/2+noise*eps;
            case 15 %Two Parabolas, dim=1, n=100
                y=( xA.^2  + unifrnd(0,1,n,1))/2.*(randsample(2, n, true)-1.5)*2+noise*eps;
            case 16 %Circle, dim=1, n=100
                x=sin(xA*pi)+mvnrnd(0,1,n)/8;
                y=cos(xA*pi)+mvnrnd(0,1,n)/8+noise*eps;
            case 17 %Independent clouds, dim=1, n=100
                x=mvnrnd(zeros(n,d),eye(d),n)*A/3+(randsample(2, n, true)-1.5)*2;
                y=mvnrnd(zeros(n,d),eye(d),n)*A/3+(randsample(2, n, true)-1.5)*2+noise*eps;
            case 18 %Joint Normal, dim=5, n=200
                rho=0.1;
                cov1=[eye(d) rho*ones(d)];
                cov2=[rho*ones(d) eye(d)];
                %                     cov1=[eye(d) rho*eye(d)];
                %                     cov2=[rho*eye(d) eye(d)];
                cov=[cov1' cov2'];
                x=mvnrnd(zeros(n,2*d),cov,n);
                y=x(:,d+1:2*d)+noise*repmat(eps,1,d);
                x=x(:,1:d);
            case 19 %Log(X^2), dim=5, n=100
                x=mvnrnd(zeros(n, d),eye(d));
                y=log(x.^2)+noise*repmat(eps,1,d);
            case 20 %Multiplicative Noise, dim=5, n=100
                x=mvnrnd(zeros(n, d),eye(d));
                y=mvnrnd(zeros(n, d),eye(d));
                y=x.*y+noise*repmat(eps,1,d);
        end
    end
    
    % Iterate through all sample size choices
    for i=1:lim
        % Form the distance matrix
        nn=numRange(i);
        C=squareform(pdist(x(1:nn,:)));
        P=squareform(pdist(y(1:nn,:)));
        
        % Permute the second dataset for rep2 times, and calculate the 4 correlation statistics
        for r2=1:rep2
            per=randperm(nn);
            Pa=P(per,per);
            dCor1(r2)=distCorr(C,Pa);
            dCor2(r2,1:nn-1)=rankDCorr(C,Pa);
            dCor3(r2)=mDistCorr(C,Pa);
            dCor4(r2)=HHG(C,Pa);
        end
        
        % Calculate the correlation statistics for the original data, and
        % derive the p-value for the permutation test.
        cut1=distCorr(C,P);
        power1(i,r1)=mean(dCor1<cut1);
        cut2=rankDCorr(C,P);
        power2(i,1:nn-1,r1)=mean(dCor2(:,1:nn-1)<repmat(cut2,1,rep2)',1);
        cut3=mDistCorr(C,P);
        power3(i,r1)=mean(dCor3<cut3);
        cut4=HHG(C,P);
        power4(i,r1)=mean(dCor4<cut4);
    end
end

% Output the mean p-value and standard deviation. Std is not meaningful
% when rep1=1.
std1=std(power1,0,2);
std2=std(power2,0,3);
std3=std(power3,0,2);
std4=std(power4,0,2);
power1=1-mean(power1,2);
power2=1-mean(power2,3);
power3=1-mean(power3,2);
power4=1-mean(power4,2);

% Save the results
if size(type,1)>1
    type=' Brain CxP';
end
filename=strcat('CorrPermTestType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
save(filename,'power1','power2','power3','power4','std1','std2','std3','std4','dCor1','dCor2','dCor3','dCor4','cut1','cut2','cut3','cut4','type','n','numRange','dim','noise','lim','rep1','rep2');

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
            titlechar=' Log(X+1)';
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
    
    % Plot figure
    xaxis=1:lim;
    % Plot the best power among all neighborhood choice of rdCorr
    plot(xaxis, power1,'b*-',xaxis,min(power2,[],2),'ro-',xaxis,power3,'g+-',xaxis,power4,'cs-','LineWidth',2);
    % Plot the power at the largest neighborhood size of rdCorr
    %plot(xaxis, power1,'b*-',xaxis,power2(:,end),'ro-',xaxis,power3,'g+-',xaxis,power4,'cs-','LineWidth',2);
    legend('Distance Correlation','Ranking Distance Correlation','Modified Distance Correlation','HHG','Location','NorthEast');
    
    % Figure title/labels
    xlabel('Sample Size');
    ylabel('Empirical Testing Power');
    xlim([1 lim]);
    ylim([0 1]);
    ax=gca;
    ax.XTickLabel=numRange;
    titleStr = strcat('Permutation Test for ', titlechar,' at n=', num2str(n), ' and m= ',num2str(dim));
    title(titleStr);
    saveas(gcf,filename,'jpeg');
end
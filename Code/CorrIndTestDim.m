function [power1, power2, power3, power4]=CorrIndTestDim(type,n,dim,lim,rep1,rep2)
% Author: Cencheng Shen
% Independence Tests for identifying correlation, returning empirical testing power with respect to increasing dimension.
% The output is the powers of rankDCorr, dCorr, modified dCorr, HHG.
%
% Type 1-10: for Y=f(XA), similar to Tibs comment 2008
% n=100; dim=1000; lim=10; rep1=1; rep2=1000;
% CorrIndTestDim(1,n,dim,lim,rep1,rep2);
%
% Type 11-17: seven correlation types in Wikipedia correlation figure /
% Table 3 of HHG, which uses dim=1.
% n=100; dim=50; lim=10; rep1=1; rep2=1000;
% CorrIndTestDim(11,n,dim,lim,rep1,rep2);
%
% Type 18-20: 3 correlation types in Szeley 2007, which uses dim=5.
% n=100; dim=10; lim=10; rep1=1; rep2=1000;
% CorrIndTestDim(18,n,dim,lim,rep1,rep2);

% Parameters:
% In the input, lim specifies the number of intervals in the dimension,
% Rep specifies the number of MC-replicates.
display=0; %change to 1 for picture auto-saving
noise=0; %noise level, if any
K=n;
alpha=0.05; %type 1 error level
option1=1; option2=1; option3=1; option4=1; %Control whether to calculate the respective correlation statistic or not.

if lim==0
    dimRange=dim;
    lim=1;
else
    dimRange=ceil(dim/lim):ceil(dim/lim):dim; %do dimension choice at the interval of ceil(dim/lim).
end

% High-dimensional decay
AA=ones(dim,1);
%AA=AA./(ceil(dim*rand(dim,1))); %random decay
for d=1:dim
    AA(d)=AA(d)/d; %fixed decay
end

% Output
power1=zeros(lim,K,rep1);power2=zeros(lim,K,rep1);power3=zeros(lim,K,rep1);power4=zeros(lim,rep1);%powers for rankdCorr, dCorr, modified dCorr, HHG
dCor1A=zeros(rep2,K,lim);dCor2A=zeros(rep2,K,lim);dCor3A=zeros(rep2,K,lim);dCor4A=zeros(rep2,lim);
dCor1N=zeros(rep2,K,lim);dCor2N=zeros(rep2,K,lim);dCor3N=zeros(rep2,K,lim);dCor4N=zeros(rep2,lim);

%rep1 MC replicates
for r1=1:rep1
    % First generate alternative distribution, i.e., X and Y are independent
    for r2=1:rep2
        % Iterate through all dimension range
        for i=1:lim
            d=dimRange(i);
            A=AA(1:d);
            x=unifrnd(-1,1,n,d);
            %x=tan((2*x-1)*pi/2); %Cauchy
            %x=mvnrnd(zeros(n, d),eye(d));
            xA=x*A;
            
            if noise ~=0
                eps=mvnrnd(0,1,n); %Gaussian noise
            else
                eps=zeros(n,1);
            end
            x=unifrnd(-1,1,n,d);
            
            switch type %Various simulations and their example parameters
                case 1 %Linear, dim=100, n=100
                    y=xA+2*noise*eps;
                case 2 %Quadratic, dim=100, n=200
                    y=4*(xA-0.5).^2+10*noise*eps;
                case 3 %Cubic, dim=100, n=200
                    y=128*(xA-1/3).^3+48*(xA-1/3).^2-12*(xA-1/3)+500*noise*eps;
                case 4 %Sine 1/2, dim=1, n=100
                    y=sin(4*pi*xA/max(xA))+2*noise*eps;
                case 5 %Sine 1/8, dim=1, n=200
                    y=sin(16*pi*xA/max(xA))+noise*eps;
                    %y=(xA-1/3).^10+1000*noise*eps;
                case 6 %X^0.25, dim=1, n=50
                    y=abs(xA).^0.25+noise/2*eps;
                case 7 %Circle, dim=10, n=100
                    y=(binornd(1,0.5,n,1)*2-1).*abs(1-(xA/max(xA)*2-1).^2).^0.5+noise*eps;
                case 8 %Step Function, dim=100, n=100
                    y=(xA>mean(xA))+1*noise*eps;
                case 9 %Exponential, dim=100, n=100
                    y=exp(xA)+2*noise*eps;
                case 10 %Uncorrelated X & Y, dim=1, n=40
                    x=binornd(1,0.5,n,d);
                    y=(binornd(1,0.5,n,d)*2-1)*A;
                    y=x*A.*y+1*noise*eps;
                    x=binornd(1,0.5,n,d);
                case 11 %W Shape, dim=1, n=100
                    y=4*( ( xA.^2 - 1/2 ).^2 + unifrnd(0,1,n,d)*A/500 )+noise*eps;
                    x=x+unifrnd(0,1,n,d)/3;
                case 12 %Square, dim=1, n=200
                    u=unifrnd(-1,1,n,d)*A;
                    v=unifrnd(-1,1,n,d)*A;
                    theta=-pi/8;
                    tmp=[cos(theta) -sin(theta); sin(theta) cos(theta)];
                    uv=[u v] * tmp;
                    y=uv(:,2)+noise*eps;
                    u=unifrnd(-1,1,n,d)*A;
                    v=unifrnd(-1,1,n,d)*A;
                    uv=[u v] * tmp;
                    x=uv(:,1);
                case 13 %Diamond, dim=1, n=200
                    u=unifrnd(-1,1,n,d)*A;
                    v=unifrnd(-1,1,n,d)*A;
                    theta=-pi/4;
                    tmp=[cos(theta) -sin(theta); sin(theta) cos(theta)];
                    uv=[u v] * tmp;
                    y=uv(:,2)+noise*eps;
                    u=unifrnd(-1,1,n,d)*A;
                    v=unifrnd(-1,1,n,d)*A;
                    uv=[u v] * tmp;
                    x=uv(:,1);
                case 14 %Parabola, dim=1, n=100
                    y=( xA.^2  + unifrnd(0,1,n,1))/2+noise*eps;
                case 15 %Two Parabolas, dim=1, n=100
                    y=( xA.^2  + unifrnd(0,1,n,1))/2.*(randsample(2, n, true)-1.5)*2+noise*eps;
                case 16 %Circle, dim=1, n=100
                    y=cos(xA*pi)+mvnrnd(0,1,n)/8+noise*eps;
                    x=sin(x*A*pi)+mvnrnd(0,1,n)/8;
                case 17 %Independent clouds, dim=1, n=100
                    x=mvnrnd(zeros(n,d),eye(d),n)*A/3+(randsample(2, n, true)-1.5)*2;
                    y=mvnrnd(zeros(n,d),eye(d),n)*A/3+(randsample(2, n, true)-1.5)*2+noise*eps;
                case 18 %Joint Normal, dim=5, n=200
                    rho=1/(dim*2);
                    cov1=[eye(d) rho*ones(d)];
                    cov2=[rho*ones(d) eye(d)];
                    %                 cov1=[eye(d) rho*eye(d)];
                    %                 cov2=[rho*eye(d) eye(d)];
                    cov=[cov1' cov2'];
                    x=mvnrnd(zeros(n,2*d),cov,n);
                    y=x(:,d+1:2*d)+2*noise*repmat(eps,1,d);
                    x=mvnrnd(zeros(n,2*d),cov,n);
                    x=x(:,1:d);
                case 19 %Log(X^2), dim=5, n=100
                    x=mvnrnd(zeros(n, d),eye(d));
                    y=log(x.^2)+10*noise*repmat(eps,1,d);
                    x=mvnrnd(zeros(n, d),eye(d));
                case 20 %Multiplicative Noise, dim=5, n=100
                    x=mvnrnd(zeros(n, d),eye(d));
                    y=mvnrnd(zeros(n, 1),eye(1));
                    y=x.*repmat(y,1,d)+2*noise*repmat(eps,1,d);
                    x=mvnrnd(zeros(n, d),eye(d));
            end
            
            % Form the distance matrix and calculate the 4 correlation statistics under the alternative
            C=squareform(pdist(x));
            P=squareform(pdist(y));
            disRank=[disToRanks(C) disToRanks(P)];
            if option1~=0
                dCor1A(r2,:,i)=rankDCorr(C,P,disRank,1);
            end
            if option2~=0
                dCor2A(r2,:,i)=rankDCorr(C,P,disRank,2);
            end
            if option3~=0
                dCor3A(r2,:,i)=rankDCorr(C,P,disRank,3);
            end
            if option4~=0
                dCor4A(r2,i)=HHG(C,P);
            end
        end
    end
    
    % Then generate null distribution, i.e., X and Y are correlated by certain
    % function type.
    for r2=1:rep2
        for i=1:lim
            d=dimRange(i);
            A=AA(1:d);
            x=unifrnd(-1,1,n,d);
            %x=tan((2*x-1)*pi/2);%Cauchy
            %x=mvnrnd(zeros(n, dim),eye(dim));
            xA=x*A;
            
            if noise ~=0
                eps=mvnrnd(0,1,n);
            else
                eps=zeros(n,1);
            end
            
            switch type %Various simulations and their example parameters
                case 1 %Linear, dim=100, n=100
                    y=xA+2*noise*eps;
                case 2 %Quadratic, dim=100, n=200
                    y=4*(xA-0.5).^2+10*noise*eps;
                case 3 %Cubic, dim=100, n=200
                    y=128*(xA-1/3).^3+48*(xA-1/3).^2-12*(xA-1/3)+500*noise*eps;
                case 4 %Sine 1/2, dim=1, n=100
                    y=sin(4*pi*xA/max(xA))+2*noise*eps;
                case 5 %Sine 1/8, dim=1, n=200
                    y=sin(16*pi*xA/max(xA))+noise*eps;
                    %y=(xA-1/3).^10+1000*noise*eps;
                case 6 %X^0.25, dim=1, n=40
                    y=abs(xA).^0.25+noise/2*eps;
                case 7 %Circle, dim=10, n=100
                    y=(binornd(1,0.5,n,1)*2-1).*abs(1-(xA/max(xA)*2-1).^2).^0.5+noise*eps;
                case 8 %Step Function, dim=100, n=100
                    y=(xA>mean(xA))+1*noise*eps;
                case 9 %Exponential, dim=100, n=100
                    y=exp(xA)+2*noise*eps;
                case 10 %Uncorrelated X & Y, dim=1, n=40
                    x=binornd(1,0.5,n,d);
                    y=(binornd(1,0.5,n,d)*2-1)*A;
                    y=x*A.*y+1*noise*eps;
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
                    y=log(x.^2)+10*noise*repmat(eps,1,d);
                case 20 %Multiplicative Noise, dim=5, n=100
                    x=mvnrnd(zeros(n, d),eye(d));
                    y=mvnrnd(zeros(n, 1),eye(1));
                    y=x.*repmat(y,1,d)+2*noise*repmat(eps,1,d);
            end
            
            % Form the distance matrix and calculate the 4 correlation
            % statistics under the null
            C=squareform(pdist(x));
            P=squareform(pdist(y));
            disRank=[disToRanks(C) disToRanks(P)];
            if option1~=0
                dCor1N(r2,:,i)=rankDCorr(C,P,disRank,1);
            end
            if option2~=0
                dCor2N(r2,:,i)=rankDCorr(C,P,disRank,2);
            end
            if option3~=0
                dCor3N(r2,:,i)=rankDCorr(C,P,disRank,3);
            end
            if option4~=0
                dCor4N(r2,i)=HHG(C,P);
            end
        end
    end
    
    % Based on the emprical test statistics under the null and the alternative,
    % calculate the emprical power at type 1 error level alpha=0.05
    for i=1:lim
        for k=1:K
            dCorT=sort(dCor1A(:,k,i),'descend');
            cut1=dCorT(ceil(rep2*alpha));
            power1(i,k,r1)=mean(dCor1N(:,k,i)>cut1);
            
            dCorT=sort(dCor2A(:,k,i),'descend');
            cut2=dCorT(ceil(rep2*alpha));
            power2(i,k,r1)=mean(dCor2N(:,k,i)>cut2);
            
            dCorT=sort(dCor3A(:,k,i),'descend');
            cut3=dCorT(ceil(rep2*alpha));
            power3(i,k,r1)=mean(dCor3N(:,k,i)>cut3);
        end
        dCorT=sort(dCor4A(:,i),'descend');
        cut4=dCorT(ceil(rep2*alpha));
        power4(i,r1)=mean(dCor4N(:,i)>cut4);
    end
end
std1=std(power1,0,3);
std2=std(power2,0,3);
std3=std(power3,0,3);
std4=std(power4,0,2);
power1=mean(power1,3);
power2=mean(power2,3);
power3=mean(power3,3);
power4=mean(power4,2);

% Save the results
filename=strcat('CorrIndTestDimType',num2str(type),'N',num2str(n),'Dim');
%save(filename,'power1','power2','power3','power4','type','n','dimRange','std1','std2','std3','std4','rep1','rep2','lim','dim','dCor1N','dCor2N','dCor3N','dCor4N','dCor1A','dCor2A','dCor3A','dCor4A');
save(filename,'power1','power2','power3','power4','type','n','dimRange','std1','std2','std3','std4','rep1','rep2','lim','dim');

% Display and save picture. By default do not display.
if display~=0
    filename=strcat('CorrIndTestDimType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
    titlechar=' Data';
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
    
    % Plot figure
    % Plot the best power among all neighborhood choice
    %plot(dimRange,max(power1,[],2),'ro-',dimRange, max(power2,[],2),'bx-',dimRange,max(power3,[],2),'g^-',dimRange,power4,'c.-','LineWidth',2);
    % Plot the power at the largest neighborhood size of dCorr
    plot(dimRange, max(power1,[],2),'ro-',dimRange,power2(:,end),'bx-',dimRange,power3(:,end),'g^-',dimRange,power4,'c.-','LineWidth',2);
    legend('Rank Distance Correlation','Distance Correlation','Modified Distance Correlation','HHG','Location','SouthEast');
    
    % Figure title/labels
    xlabel('Dimension');
    ylabel('Empirical Testing Power');
    xlim([dimRange(1) dimRange(end)]);
    ylim([0 1]);
    titleStr = strcat('Independence Test for ', titlechar,' at n=', num2str(n), ' up To m= ',num2str(dim));
    title(titleStr);
    saveas(gcf,filename,'jpeg');
end
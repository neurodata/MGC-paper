function [power1, power2, power3, power4]=CorrIndTestDim(type,n,dim,lim,rep)
% Author: Cencheng Shen
% Independence Tests for identifying correlation, returning empirical testing power with respect to increasing dimension.
% The output is the powers of dCorr, rankDCorr, modified dCorr, HHG.
%
% Type 1-10: for Y=f(XA), similar to Tibs comment 2008
% n=100; dim=500; lim=10; rep=1000;
% CorrIndTestDim(1,n,dim,lim,rep);
%
% Type 11-17: seven correlation types in Wikipedia correlation figure /
% Table 3 of HHG, which uses dim=1.
% n=100; dim=50; lim=10; rep=1000;
% CorrIndTestDim(11,n,dim,lim,rep);
%
% Type 18-20: 3 correlation types in Szeley 2007, which uses dim=5.
% n=100; dim=10; lim=10; rep=1000;
% CorrIndTestDim(18,n,dim,lim,rep);

% Parameters:
% In the input, lim specifies the number of intervals in the dimension,
% Rep specifies the number of MC-replicates.
display=0; %change to 1 for picture auto-saving
noise=0; %noise level, if any
K=n-1;
alpha=0.05; %type 1 error level

if lim==0
    dimRange=dim;
    lim=1;
else
    dimRange=0:ceil(dim/lim):dim; %do dimension at the interval of lim from 1 to dim.
    dimRange(1)=1;
    lim=length(dimRange);
end

% High-dimensional decay
AA=ones(dim,1);
%AA=AA./(ceil(dim*rand(dim,1))); %random decay
for d=1:dim
    AA(d)=AA(d)/d; %fixed decay
end

% Output
power1=zeros(lim,1); power2=zeros(lim,K);power3=zeros(lim,1);power4=zeros(lim,1);%powers for dCorr, rankdCorr, modified dCorr, HHG
dCor1A=zeros(rep,lim);dCor2A=zeros(rep,K,lim);dCor3A=zeros(rep,lim);dCor4A=zeros(rep,lim);
dCor1N=zeros(rep,lim);dCor2N=zeros(rep,K,lim);dCor3N=zeros(rep,lim);dCor4N=zeros(rep,lim);

% First generate alternative distribution, i.e., X and Y are independent
for r=1:rep
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
            eps=0;
        end
        x=unifrnd(-1,1,n,d);
        
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
                rho=0.1;
                cov1=[eye(d) rho*ones(d)];
                cov2=[rho*ones(d) eye(d)];
                %                 cov1=[eye(d) rho*eye(d)];
                %                 cov2=[rho*eye(d) eye(d)];
                cov=[cov1' cov2'];
                x=mvnrnd(zeros(n,2*d),cov,n);
                y=x(:,d+1:2*d)+noise*repmat(eps,1,d);
                x=mvnrnd(zeros(n,2*d),cov,n);
                x=x(:,1:d);
            case 19 %Log(X^2), dim=5, n=100
                x=mvnrnd(zeros(n, d),eye(d));
                y=log(x.^2)+noise*repmat(eps,1,d);
                x=mvnrnd(zeros(n, d),eye(d));
            case 20 %Multiplicative Noise, dim=5, n=100
                x=mvnrnd(zeros(n, d),eye(d));
                y=mvnrnd(zeros(n, d),eye(d));
                y=x.*y+noise*repmat(eps,1,d);
                x=mvnrnd(zeros(n, d),eye(d));
        end
        
        % Form the distance matrix and calculate the 4 correlation statistics under the alternative
        C=squareform(pdist(x));
        P=squareform(pdist(y));
        dCor1A(r,i)=distCorr(C,P);
        dCor2A(r,:,i)=rankDCorr(C,P);
        dCor3A(r,i)=mDistCorr(C,P);
        dCor4A(r,i)=HHG(C,P);
    end
end

% Then generate null distribution, i.e., X and Y are correlated by certain
% function type.
for r=1:rep
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
            eps=0;
        end
        
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
        
        % Form the distance matrix and calculate the 4 correlation
        % statistics under the null
        C=squareform(pdist(x));
        P=squareform(pdist(y));
        dCor1N(r,i)=distCorr(C,P);
        dCor2N(r,:,i)=rankDCorr(C,P);
        dCor3N(r,i)=mDistCorr(C,P);
        dCor4N(r,i)=HHG(C,P);
    end
end

% Based on the emprical test statistics under the null and the alternative,
% calculate the emprical power at type 1 error level alpha=0.05
for i=1:lim
    dCorT=sort(dCor1A(:,i),'descend');
    cut1=dCorT(ceil(rep*alpha));
    dCorT=sort(dCor3A(:,i),'descend');
    cut3=dCorT(ceil(rep*alpha));
    dCorT=sort(dCor4A(:,i),'descend');
    cut4=dCorT(ceil(rep*alpha));
    power1(i)=mean(dCor1N(:,i)>cut1);
    power3(i)=mean(dCor3N(:,i)>cut3);
    power4(i)=mean(dCor4N(:,i)>cut4);
    for k=1:K
        dCorT=sort(dCor2A(:,k,i),'descend');
        cut2=dCorT(ceil(rep*alpha));
        power2(i,k)=mean(dCor2N(:,k,i)>cut2);
    end
end

% Save the results
filename=strcat('CorrIndTestDimType',num2str(type),'N',num2str(n),'Dim',num2str(dim));
save(filename,'power1','power2','power3','power4','type','n','dimRange','rep','lim','dim','dCor1N','dCor2N','dCor3N','dCor4N','dCor1A','dCor2A','dCor3A','dCor4A');

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
    plot(xaxis, power1,'b*-',xaxis,max(power2(:,:),[],2),'ro-',xaxis,power3,'g+-',xaxis,power4,'cs-','LineWidth',2);
    % Plot the power at the largest neighborhood size of rdCorr
    %plot(xaxis, power1,'b*-',xaxis,power2(:,end),'ro-',xaxis,power3,'g+-',xaxis,power4,'cs-','LineWidth',2);
    legend('Distance Correlation','Ranking Distance Correlation','Modified Distance Correlation','HHG','Location','NorthEast');
    
    % Figure title/labels
    xlabel('Dimension');
    ylabel('Empirical Testing Power');
    xlim([1 lim]);
    ylim([0 1]);
    ax=gca;
    %     ax.XTick=1:(lim-1)/10:lim;
    %     dimRange2=0:dim/10:dim;
    %     dimRange2(1)=1;
    ax.XTickLabel=dimRange;
    titleStr = strcat('Independence Test for ', titlechar,' at n=', num2str(n), ' up To m= ',num2str(dim));
    title(titleStr);
    saveas(gcf,filename,'jpeg');
end
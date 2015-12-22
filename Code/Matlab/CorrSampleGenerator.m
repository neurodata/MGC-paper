function [x, y]=CorrSampleGenerator(type,n,dim,optionAlt, noise)
% Author: Cencheng Shen
% Output sample data x and y for testing independence based on given type.
%
% Parameters:
% type specifies the type of distribution,
% n is the sample size, dim is the dimension
% optionAlt specifies the null or the alternative hypothesis, by default 1
% noise specifies the noise level, by default 0.
if nargin<4
    optionAlt=1; %by default generate dependent samples
end
if nargin<5
    noise=0; %default noise level
end

if noise~=0
    eps=mvnrnd(0,1,n); %Gaussian noise added to Y
else
    eps=zeros(n,1);
end

% High-dimensional decay
A=ones(dim,1);
%A=A./(ceil(dim*rand(dim,1))); %random decay
for d=1:dim
    A(d)=A(d)/d; %fixed decay
end

% Generate x by uniform distribution first, which is the distribution used by half
% of the distributions; store the weighted summation in xA.
x=unifrnd(-1,1,n,d);
%x=random('norm',0,1,n,d);
xA=x*A;
% Generate x independently by uniform if the null hypothesis is true, i.e., x is independent of y.
if optionAlt==0
    x=unifrnd(-1,1,n,d);
    %x=random('norm',0,1,n,d);
end

switch type % In total 20 types of dependency
    case 1 %Linear
        y=xA+1*noise*eps;
    case 2 %Quadratic
        y=(xA).^2+0.5*noise*eps;
    case 3 %Cubic
        y=128*(xA-1/3).^3+48*(xA-1/3).^2-12*(xA-1/3)+80*noise*eps;
    case 4 %Sine 1/2
        y=sin(4*pi*xA/max(xA))+1*noise*eps;
    case 5 %Sine 1/8
        y=sin(16*pi*xA/max(xA))+0.5*noise*eps;
    case 6 %X^0.25
        y=abs(xA).^0.25+noise/5*eps;
    case 7 %Circle
        x=unifrnd(0,1,n,d);
        xA=x*A;
        y=(binornd(1,0.5,n,1)*2-1).*abs(1-(xA/max(xA)*2-1).^2).^0.5+0.2*noise*eps;
        if optionAlt==0
            x=unifrnd(0,1,n,d);
        end
    case 8 %Step Function
        y=(xA>mean(xA))+1*noise*eps;
    case 9 %Exponential
        y=exp(xA)+1.5*noise*eps;
    case 10 %Uncorrelated Binomial
        x=binornd(1,0.5,n,d);
        y=(binornd(1,0.5,n,d)*2-1)*A;
        y=x*A.*y+0.6*noise*eps;
        if optionAlt==0
            x=binornd(1,0.5,n,d);
        end
    case 11 %W Shape
        y=4*( ( xA.^2 - 1/2 ).^2 + unifrnd(0,1,n,d)*A/500 )+0.5*noise*eps;
    case 12 %Square Shape
        u=unifrnd(-1,1,n,d)*A;
        v=unifrnd(-1,1,n,d)*A;
        theta=-pi/8;
        tmp=[cos(theta) -sin(theta); sin(theta) cos(theta)];
        uv=[u v] * tmp;
        x=uv(:,1);
        y=uv(:,2)+noise*eps;
        if optionAlt==0
            u=unifrnd(-1,1,n,d)*A;
            v=unifrnd(-1,1,n,d)*A;
            uv=[u v] * tmp;
            x=uv(:,1);
        end
    case 13 %Diamond Shape
        u=unifrnd(-1,1,n,d)*A;
        v=unifrnd(-1,1,n,d)*A;
        theta=-pi/4;
        tmp=[cos(theta) -sin(theta); sin(theta) cos(theta)];
        uv=[u v] * tmp;
        x=uv(:,1);
        y=uv(:,2)+noise*eps;
        if optionAlt==0
            u=unifrnd(-1,1,n,d)*A;
            v=unifrnd(-1,1,n,d)*A;
            uv=[u v] * tmp;
            x=uv(:,1);
        end
    case 14 %Parabola
        y=( xA.^2  + unifrnd(0,1,n,1))/2+0.2*noise*eps;
    case 15 %Two Parabolas
        y=( xA.^2  + unifrnd(0,1,n,1)).*(binornd(1,0.5,n,1)-0.5)+noise*eps;
    case 16 %Circle 2
        z=unifrnd(-1,1,n,d);
        y=cos(z*pi)*A+mvnrnd(0,1,n)/8+0.2*noise*eps;
        if optionAlt==0
            z=unifrnd(-1,1,n,d);
        end
        x=sin(z*pi)*A+mvnrnd(0,1,n)/8;
    case 17 %Independent clouds
        x=mvnrnd(zeros(n,d),eye(d),n)*A/3+(binornd(1,0.5,n,1)-0.5)*2;
        y=mvnrnd(zeros(n,d),eye(d),n)*A/3+(binornd(1,0.5,n,1)-0.5)*2+noise*eps;
    case 18 %Joint Normal; note that the covariance matrix for dim>10 is no longer positive semi-definite
        rho=1/(d*2);
        cov1=[eye(d) rho*ones(d)];
        cov2=[rho*ones(d) eye(d)];
        cov=[cov1' cov2'];
        x=mvnrnd(zeros(n,2*d),cov,n);
        y=x(:,d+1:2*d)+2*noise*repmat(eps,1,d);
        if optionAlt==0
            x=mvnrnd(zeros(n,2*d),cov,n);
        end
        x=x(:,1:d);
    case 19 %Log(X^2)
        x=mvnrnd(zeros(n, d),eye(d));
        y=log(x.^2)+3*noise*repmat(eps,1,d);
        if optionAlt==0
            x=mvnrnd(zeros(n, d),eye(d));
        end
    case 20 %Multiplicative Noise
        x=mvnrnd(zeros(n, d),eye(d));
        y=mvnrnd(zeros(n, 1),eye(1));
        y=x.*repmat(y,1,d)+1*noise*repmat(eps,1,d);
        if optionAlt==0
            x=mvnrnd(zeros(n, d),eye(d));
        end
end
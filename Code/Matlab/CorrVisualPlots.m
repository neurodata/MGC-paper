function []=CorrVisualPlots(n,dim,noise)
% Author: Cencheng Shen
% n=1000;dim=1;noise=0;
% CorrVisualPlots(n,dim,noise)
% Used to plot figure 0 in the files
A=ones(dim,1);
%A=A./(ceil(dim*rand(dim,1))); %random decay
for d=1:dim
    A(d)=A(d)/d; %fixed decay
end

d=dim;

if noise~=0
    eps=mvnrnd(0,1,n); %Gaussian noise
else
    eps=zeros(n,1);
end

figure
s=5;
t=4;
for type=1:20
    subplot(s,t,type)
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
    
    switch type %Various simulations and their example parameters
        case 1 %Linear, dim=100, n=100
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=xA+1*noise*eps;
        case 2 %Quadratic, dim=100, n=200
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=(xA).^2+0.5*noise*eps;
        case 3 %Cubic, dim=100, n=200
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=128*(xA-1/3).^3+48*(xA-1/3).^2-12*(xA-1/3)+80*noise*eps;
        case 4 %Sine 1/2, dim=1, n=100
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=sin(4*pi*xA/max(xA))+1*noise*eps;
        case 5 %Sine 1/8, dim=1, n=200
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=sin(16*pi*xA/max(xA))+0.5*noise*eps;
            %y=(xA-1/3).^10+1000*noise*eps;
        case 6 %X^0.25, dim=1, n=40
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=abs(xA).^0.25+noise/5*eps;
        case 7 %Circle, dim=10, n=100
            x=unifrnd(0,1,n,d);
            xA=x*A;
            y=(binornd(1,0.5,n,1)*2-1).*abs(1-(xA/max(xA)*2-1).^2).^0.5+noise*eps;
        case 8 %Step Function, dim=100, n=100
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=(xA>mean(xA))+1*noise*eps;
        case 9 %Exponential, dim=100, n=100
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=exp(xA)+1.5*noise*eps;
        case 10 %Uncorrelated X & Y, dim=1, n=40
            x=binornd(1,0.5,n,d);
            y=(binornd(1,0.5,n,d)*2-1)*A;
            y=x*A.*y+0.6*noise*eps;
        case 11 %W Shape, dim=1, n=100
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=4*( ( xA.^2 - 1/2 ).^2 + unifrnd(0,1,n,d)*A/500 )+0.5*noise*eps;
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
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=( xA.^2  + unifrnd(0,1,n,1))/2+0.2*noise*eps;
        case 15 %Two Parabolas, dim=1, n=100
            x=unifrnd(-1,1,n,d);
            xA=x*A;
            y=( xA.^2  + unifrnd(0,1,n,1)).*(binornd(1,0.5,n,1)-0.5)+noise*eps;
        case 16 %Circle 2, dim=1, n=100
%             x=unifrnd(-1,1,n,d);
%             xA=x*A;
%             x=sin(xA*pi)+mvnrnd(0,1,n)/8;
%             y=cos(xA*pi)+mvnrnd(0,1,n)/8+0.2*noise*eps;
            z=unifrnd(-1,1,n,d);
            %xA=x*A;
            x=sin(z*pi)*A+mvnrnd(0,1,n)/8;
            y=cos(z*pi)*A+mvnrnd(0,1,n)/8+0.2*noise*eps;
        case 17 %Independent clouds, dim=1, n=100
            x=mvnrnd(zeros(n,d),eye(d),n)*A/3+(binornd(1,0.5,n,1)-0.5)*2;
            y=mvnrnd(zeros(n,d),eye(d),n)*A/3+(binornd(1,0.5,n,1)-0.5)*2+noise*eps;
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
            y=log(x.^2)+3*noise*repmat(eps,1,d);
        case 20 %Multiplicative Noise, dim=5, n=100
            x=mvnrnd(zeros(n, d),eye(d));
            y=mvnrnd(zeros(n, 1),eye(1));
            y=x.*repmat(y,1,d)+1*noise*repmat(eps,1,d);
    end
    plot(x,y,'.')
    title(titlechar);
end
suptitle('Visualization for 20 Simulated Dependencies')

function [score,ind,score2]=MGCGeometry(A,B,rep)

% if nargin<2
%     localCorr=A;
% else
%     localCorr=MGCLocalCorr(A,B);
% end
if nargin<3
    rep=10;
end
if size(A,1)==size(A,2) && sum(diag(A).^2)==0
    x=cmdscale(A,1);
%     disp('The first data is not a Euclidean distance matrix; transformed to distance matrix instead.')
else
    x=A;
end
localCorr=MGCLocalCorr(A,B);
%
%
% ind2=zeros(K*3,1);
total=20;
n=size(localCorr,1);
% localCorr2=zeros(n,n);
    score=zeros(total,1);
%     x=cmdscale(A,1);
for k=1:rep
    for j=1:total
                                    [y]=CorrSampleGeneratorX(j,x,0);
%         for k=1:K
%         [x,y]=CorrSampleGenerator(j,n,1,1, 0.001);
        [~, tmp, ~]=MGCSampleStat(x,y);
        %     [~, ~,~,localCorr2]=MGCPermutationTest(x,y,rep);
        %     localCorr2=double(localCorr2<0.05);
        %             localCorr2=localCorr2-mean(mean(localCorr2));
        %     var2=norm(localCorr2,'fro');
        %     localCorr2=reshape(localCorr2,size(localCorr2,1)*size(localCorr2,2),1);
        %     pcorr(j)=sum(sum(localCorr.*localCorr2))/var1/var2;
%         localCorr2=localCorr2+tmp/K;
%         end
        score(j)=score(j)+DCorr(localCorr,tmp)/rep;
        %     pcorr(j)=MGCSampleStat(localCorr,localCorr2);
        %     pcorr(j)=corr(localCorr,localCorr2,'Spearman');
    end
end
    
[~,ind]=sort(score,'descend');


if size(A,1)==size(A,2) && sum(diag(A).^2)==0
    x=cmdscale(A,1);
%     disp('The first data is not a Euclidean distance matrix; transformed to distance matrix instead.')
else
    x=A;
end
if size(B,1)==size(B,2) && sum(diag(B).^2)==0
    y=cmdscale(B,1);
%     disp('The first data is not a Euclidean distance matrix; transformed to distance matrix instead.')
else
    y=B;
end
score2=zeros(total,1);
    for j=1:total
        [z]=CorrSampleGeneratorX(j,x,0);
%         z=x;
        z=z*corr(z,y)/std(z);
        z=z+mean(y)-mean(z);
%         for k=1:K
%         [x,y]=CorrSampleGenerator(j,n,1,1, 0.001);
%         [~, tmp, ~]=MGCSampleStat(x,y);
        %     [~, ~,~,localCorr2]=MGCPermutationTest(x,y,rep);
        %     localCorr2=double(localCorr2<0.05);
        %             localCorr2=localCorr2-mean(mean(localCorr2));
        %     var2=norm(localCorr2,'fro');
        %     localCorr2=reshape(localCorr2,size(localCorr2,1)*size(localCorr2,2),1);
        %     pcorr(j)=sum(sum(localCorr.*localCorr2))/var1/var2;
%         localCorr2=localCorr2+tmp/K;
%         end
        score2(j)=norm(z-y,'fro');
        %     pcorr(j)=MGCSampleStat(localCorr,localCorr2);
        %     pcorr(j)=corr(localCorr,localCorr2,'Spearman');
    end
    [~,score2]=sort(score2,'descend');
%     ind2(k)=ind(1);
%     ind2(K+k)=ind(2);
%     ind2(2*K+k)=ind(3);
% [~,ind2]=hist(ind2,unique(ind2));
% ind2=find(ind2==t);
% % mode(ind2)
% %         if (mode(ind2)==t)
% if (ind2<=3)
%     power(t)=power(t)+1/rep;
% end
% end
% end

function [y]=CorrSampleGeneratorX(type,x,noise)
%4,5,8,16,17
[n,dim]=size(x);

eps=mvnrnd(0,1,n); % Gaussian noise added to Y

% High-dimensional decay
A=ones(dim,1);
%A=A./(ceil(dim*rand(dim,1))); %random decay
for d=1:dim
    A(d)=A(d)/d; %fixed decay
end
d=dim;

% % Generate x by uniform distribution first, which is the default distribution used by many types; store the weighted summation in xA.
% x=unifrnd(-1,1,n,d);
xA=x*A;
xA=xA/max(max(xA));
% % Generate x independently by uniform if the null hypothesis is true, i.e., x is independent of y.
% if dependent==0
%     x=unifrnd(-1,1,n,d);
% end

switch type % In total 20 types of dependency + the type 0 outlier model
    case 1 %Linear
        y=xA+1*noise*eps;
    case 2 %Exponential
%         x=unifrnd(0,3,n,d);
        y=exp(xA)+10*noise*eps;
%         if dependent==0
%             x=unifrnd(0,3,n,d);
%         end
    case 3 %Cubic
        y=128*(xA-1/3).^3+48*(xA-1/3).^2-12*(xA-1/3)+80*noise*eps;
    case 4 %Joint Normal; note that dim should be no more than 10 as the covariance matrix for dim>10 is no longer positive semi-definite
        rho=1/(d*2);
        cov1=[eye(d) rho*ones(d)];
        cov2=[rho*ones(d) eye(d)];
        covT=[cov1' cov2'];
        x=mvnrnd(zeros(n,2*d),covT,n);
        y=x(:,d+1:2*d)+0.5*noise*repmat(eps,1,d);
%         if dependent==0
%             x=mvnrnd(zeros(n,2*d),covT,n);
%         end
        x=x(:,1:d);
    case 5 %Step Function
        if dim>1
            noise=1;
        end
        y=(xA>0)+1*noise*eps;
    case 6 %Quadratic
        y=(xA).^2+0.5*noise*eps;
    case 7 %W Shape
        y=4*( ( xA.^2 - 1/2 ).^2 + unifrnd(0,1,n,d)*A/500 )+0.5*noise*eps;
    case 9 %Uncorrelated Binomial
        if d>1
            noise=1;
        end
%         x=binornd(1,0.5,n,d)+0.5*noise*mvnrnd(zeros(n,d),eye(d),n);
        y=(binornd(1,0.5,n,1)*2-1);
        y=x*A.*y+0.5*noise*eps;
%         if dependent==0
%             x=binornd(1,0.5,n,d)+0.5*noise*mvnrnd(zeros(n,d),eye(d),n);
%         end
    case 10 %Log(X^2)
%         x=mvnrnd(zeros(n, d),eye(d));
        y=log(x.^2)+3*noise*repmat(eps,1,d);
%         if dependent==0
%             x=mvnrnd(zeros(n, d),eye(d));
%         end
    case 11 %Fourth root
        y=abs(xA).^(0.25)+noise/4*eps;
    case {8,16,17} %Circle & Ecllipse & Spiral
        if d>1
            noise=1;
        end
        cc=0.4;
        if type==16
            rx=ones(n,d);
        end
        if type==17
            rx=5*ones(n,d);
        end

        if type==8
            rx=unifrnd(0,5,n,1);
            ry=rx;
            rx=repmat(rx,1,d);
            z=rx;
        else
            z=unifrnd(-1,1,n,d);
            ry=ones(n,1);
        end
        x(:,1)=cos(z(:,1)*pi);
        for i=1:d-1;
            x(:,i+1)=x(:,i).*cos(z(:,i+1)*pi);
            x(:,i)=x(:,i).*sin(z(:,i+1)*pi);
        end
        x=rx.*x;
        y=ry.*sin(z(:,1)*pi);
        if type==8
            y=y+cc*(dim)*noise*mvnrnd(zeros(n, 1),eye(1));
        else
            x=x+cc*noise*rx.*mvnrnd(zeros(n, d),eye(d));
        end
%         if dependent==0
%             if type==8
%                 rx=unifrnd(0,5,n,1);
%                 rx=repmat(rx,1,d);
%                 z=rx;
%             else
%                 z=unifrnd(-1,1,n,d);
%             end
%             x(:,1)=cos(z(:,1)*pi);
%             for i=1:d-1;
%                 x(:,i+1)=x(:,i).*cos(z(:,i+1)*pi);
%                 x(:,i)=x(:,i).*sin(z(:,i+1)*pi);
%             end
%             x=rx.*x;
%             if type==8
%             else
%                 x=x+cc*noise*rx.*mvnrnd(zeros(n, d),eye(d));
%             end
%         end
    case {12,13} %Sine 1/2 & 1/8
%         x=repmat(unifrnd(-1,1,n,1),1,d);
%         if noise>0 || d>1
%             x=x+0.02*(d)*mvnrnd(zeros(n,d),eye(d),n);
%         end
        if type==12
            theta=4;cc=1;
        else
            theta=16;cc=0.5;
        end
        y=sin(theta*pi*x)+cc*noise*repmat(eps,1,d);
%         if dependent==0
%             x=repmat(unifrnd(-1,1,n,1),1,d);
%             if noise>0 || d>1
%                 x=x+0.02*(d)*mvnrnd(zeros(n,d),eye(d),n);
%             end
%         end
    case {14,18} %Square & Diamond
        u=repmat(unifrnd(-1,1,n,1),1,d);
        v=repmat(unifrnd(-1,1,n,1),1,d);
        if type==14
            theta=-pi/8;
        else
            theta=-pi/4;
        end
        eps=0.05*(d)*mvnrnd(zeros(n,d),eye(d),n);
%         x=u*cos(theta)+v*sin(theta)+eps;
        y=-u*sin(theta)+v*cos(theta);
%         if dependent==0
%             u=repmat(unifrnd(-1,1,n,1),1,d);
%             v=repmat(unifrnd(-1,1,n,1),1,d);
%             eps=0.05*(d)*mvnrnd(zeros(n,d),eye(d),n);
%             x=u*cos(theta)+v*sin(theta)+eps;
%         end
    case 15 %Two Parabolas
        y=( xA.^2  + 2*noise*unifrnd(0,1,n,1)).*(binornd(1,0.5,n,1)-0.5);
    case 19 %Multiplicative Noise
%         x=mvnrnd(zeros(n, d),eye(d));
        y=mvnrnd(zeros(n, d),eye(d));
        y=x.*y;
%         if dependent==0
%             x=mvnrnd(zeros(n, d),eye(d));
%         end
    case 20 %Independent clouds
%         x=mvnrnd(zeros(n,d),eye(d),n)/3+(binornd(1,0.5,n,d)-0.5)*2;
        y=mvnrnd(zeros(n,d),eye(d),n)/3+(binornd(1,0.5,n,d)-0.5)*2;
end
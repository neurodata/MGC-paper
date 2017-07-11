function [dv,dc,dcr]=Screening(opt,rep,p,q)
% [dv,dc,dcr]=Screening(1,10,1,10)
% [dv,dc,dcr]=Screening(2,10,1,10)
if nargin<1
    opt=1;
end
if nargin<2
    rep=10;
end
if nargin<3
    p=1;
end
if nargin<4
    q=10;
end
dc=zeros(q+1,1);
dv=zeros(q+1,1);
dcr=zeros(q+1,1);

n=100;
x=unifrnd(0,10,n,p);
y=x;
D=squareform(pdist(y));
for q1=0:q
    %q=1;
        
    for r=1:rep
        if opt==2
            x2=unifrnd(0,1,n,q1*100);
            %x2=normrnd(0,0.2,n,q1);
            %C=squareform(pdist(x+q/q*x2));
            % if q>0
            xN=[x';x2']';
%             for j=1:size(xN,2);
%                 tmp=xN(:,j);
%                 xN(:,j)=tmp/norm(tmp);
%             end
            C=squareform(pdist(xN));
        else
            x2=unifrnd(0,0.1,n,p);
            %x2=normrnd(0,0.1,n,p);
            xN=x+q1/q*x2;
%             for j=1:size(xN,2);
%                 tmp=xN(:,j);
%                 xN(:,j)=tmp/norm(tmp);
%             end
            C=squareform(pdist(xN));
        end
        % else
        %     C=squareform(pdist(x));
        % end
        option1='dcorDouble';
        option2='mcor';
        
        [C,D,RX,RY]=MGCDistTransform(C,D,option2);
        
        dcov=sum(sum(C.*D))/n^2;
        dvar1=(sum(sum(C.*C'))/n^2)^0.5;
        dvar2=(sum(sum(D.*D'))/n^2)^0.5;
        dcorr=dcov/dvar1/dvar2;
        dv(q1+1)=dv(q1+1)+dvar1/rep;
        dc(q1+1)=dc(q1+1)+dcov/rep;
        dcr(q1+1)=dcr(q1+1)+dcorr/rep;
    end
end
x=0:q;
lw=3;
figure
hold on
plot(x,dc,'r+-','linewidth',lw);
plot(x,dcr,'m.-','linewidth',lw);
plot(x,dv,'k--','linewidth',lw);
hold off
set(gca,'yscale','log','fontSize',15)
if opt==2
xlabel('Number of Redundant Ambient Dimension');
set(gca,'XTick',[0,q/2,q]);
else
    xlabel('Amount of Noise');
    set(gca,'XTick',[0,q/2,q],'XTickLabel',[0,0.5,1]);
end
ylabel('Magnitude (Log Scale)');
legend('Distance Covariance','Distance Correlation','Distance Variance','location','southwest');
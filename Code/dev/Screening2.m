function [dcr1,dcr2,p1,p2]=Screening2(opt,rep,p,q)
% [dcr1,dcr2,p1,p2]=Screening2(1,100,1,20)
% [dcr1,dcr2,p1,p2]=Screening2(1,100,1,100)
% [dcr1,dcr2,p1,p2]=Screening2(1,100,1,500)
if nargin<1
    opt=1;
end
if nargin<2
    rep=100;
end
if nargin<3
    p=1;
end
if nargin<4
    q=10;
end

n=100;
x=unifrnd(0,10,n,p);
y=x;
D=squareform(pdist(y));
sz=50;
int=2;
repp=100;

dc=zeros((sz/int)+1,1);
dv=zeros((sz/int)+1,1);
dcr1=zeros((sz/int)+1,1);
dcr2=zeros((sz/int)+1,1);
p1=zeros((sz/int)+1,1);
p2=zeros((sz/int)+1,1);
alpha=0.05;

for q1=0:sz/int
    %q=1;
        
    for r=1:rep
%         if opt==2
%             x2=unifrnd(0,1,n,q1);
%             %C=squareform(pdist(x+q/q*x2));
%             % if q>0
%             xN=[x';x2']';
% %             for j=1:size(xN,2);
% %                 tmp=xN(:,j);
% %                 xN(:,j)=tmp/norm(tmp);
% %             end
%             C=squareform(pdist(xN));
%         else
            x2=unifrnd(0,q1*int,n,q);
            xN=[x';x2']';
            C=squareform(pdist(xN));
            for j=1:size(xN,2)
                tmp=xN(:,j);
                if norm(tmp)~=0
                xN(:,j)=tmp/norm(tmp);
                end
            end
            C2=squareform(pdist(xN));
%         end
        % else
        %     C=squareform(pdist(x));
        % end
        option1='mgc';
        option2='mcor';
        
%         [C1,D,RX,RY]=MGCDistTransform(C,D,option1);        
%         dcov=sum(sum(C1.*D1))/n^2;
%         dvar1=(sum(sum(C1.*C1'))/n^2)^0.5;
%         dvar2=(sum(sum(D1.*D1'))/n^2)^0.5;
%         dcorr=dcov/dvar1/dvar2;
        [pMGC,statMGC,pLocalCorr,localCorr,optimalInd]=MGCPermutationTest(C,D,repp,option1);
        dcorr=localCorr(end);
        pp1=pLocalCorr(end);        
%         dv(q1+1)=dv(q1+1)+dvar1/rep;
%         dc(q1+1)=dc(q1+1)+dcov/rep;
        dcr1(q1+1)=dcr1(q1+1)+dcorr/rep;
        p1(q1+1)=p1(q1+1)+(pp1<alpha)/rep;
        
        %         [C1,D1,RX,RY]=MGCDistTransform(C2,D,option1);
        %         dcov=sum(sum(C1.*D1))/n^2;
        % %         dvar1=(sum(sum(C1.*C1'))/n^2)^0.5;
        % %         dvar2=(sum(sum(D1.*D1'))/n^2)^0.5;
        %         dcorr=dcov/dvar1/dvar2;
        %         dcr2(q1+1)=dcr2(q1+1)+dcorr/rep;
        [pMGC,statMGC,pLocalCorr,localCorr,optimalInd]=MGCPermutationTest(C2,D,repp,option1);
        dcorr=localCorr(end);
        pp2=pLocalCorr(end);
        %         dv(q1+1)=dv(q1+1)+dvar1/rep;
        %         dc(q1+1)=dc(q1+1)+dcov/rep;
        dcr2(q1+1)=dcr2(q1+1)+dcorr/rep;
        p2(q1+1)=p2(q1+1)+(pp2<alpha)/rep;
    end
end
x=0:sz/int;
lw=3;
figure
hold on
dcr1(dcr1<=0)=0.001;
plot(x,dcr1,'m.-','linewidth',lw);
dcr1(dcr2<=0)=0.001;
plot(x,dcr2,'r+-','linewidth',lw);
% plot(x,dv,'k--','linewidth',lw);
hold off
set(gca,'yscale','log','fontSize',15)
% if opt==2
xlabel('Noise Level');
set(gca,'XTick',[0,10,20,30,40,50]/int,'XTickLabel',[0,1,2,3,4,5],'fontSize',15);
% else
%     xlabel('Amount of Noise');
%     set(gca,'XTick',[0,q/2,q],'XTickLabel',[0,0.5,1]);
% end
ylim([0,1]);
ylabel('Distance Correlation');
legend('Dcorr (Original)','Dcorr (Normalized)','location','southwest');
title(strcat('Dcorr at Noisy Dimension ', {' '},num2str(q)));

figure
hold on
plot(x,p1,'m.-','linewidth',lw);
plot(x,p2,'r+-','linewidth',lw);
% plot(x,dv,'k--','linewidth',lw);
hold off
% set(gca,'yscale','log','fontSize',15)
% if opt==2
xlabel('Noise Level');
set(gca,'XTick',[0,10,20,30,40,50]/int,'XTickLabel',[0,1,2,3,4,5],'fontSize',15);
% else
%     xlabel('Amount of Noise');
%     set(gca,'XTick',[0,q/2,q],'XTickLabel',[0,0.5,1]);
% end
ylabel('Testing Power');
ylim([0,1]);
legend('Dcorr (Original)','Dcorr (Normalized)','location','southwest');
title(strcat('Testing Powers at Noisy Dimension ', {' '},num2str(q)));

save(strcat('HDTrialN',num2str(n),'Dim',num2str(q)),'n','p','q','sz','int','rep','repp','p1','p2','dcr1','dcr2');
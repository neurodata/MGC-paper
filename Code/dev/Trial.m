load proteomics.mat

% per=(LabelIndAll==1 | LabelIndAll==2);
% Ovarian vs Normal
per=(LabelIndAll==1 | LabelIndAll==4);
D=LabelIndAll(per);
D=squareform(pdist(D));
D(D>0)=1;
rep=5000;
m=318;

% C=squareform(pdist(A(:,per)'));
% CorrPermDistTest(C,D,rep*10,'ProteomicsOvavsNormal');

% Screening Ovariance vs Normal
pMGC=zeros(m,1);pD=zeros(m,1);pM=zeros(m,1);pP=zeros(m,1);pHHG=zeros(m,1);testMGC=zeros(m,1);testD=zeros(m,1);testM=zeros(m,1);testHHG=zeros(m,1);
for i=1:m
    i
    C=squareform(pdist(A(i,per)'));   
    [pMGC(i),pD(i),pM(i),pP(i), pHHG(i),testMGC(i),testD(i),testM(i),testHHG(i)]=CorrPermDistTest(C,D,rep);
end

testMGC(:,6)=pMGC;
testMGC(:,7)=pP;
testMGC(:,8)=pD;
testMGC(:,9)=pM;
testMGC(:,10)=pHHG;
testMGC(:,2)=testD;
testMGC(:,3)=testD;
testMGC(:,4)=testM;
testMGC(:,5)=testHHG;
save('ScreeningOvavsNormal');


clear
load proteomics.mat
% Pancreatic vs Normal
per=(LabelIndAll==1 | LabelIndAll==2);
D=LabelIndAll(per);
D=squareform(pdist(D));
D(D>0)=1;
rep=5000;
m=318;

% C=squareform(pdist(A(:,per)'));
% CorrPermDistTest(C,D,rep*10,'ProteomicsPancvsNormal');
ind=181;%296 for hhg
C=squareform(pdist(A(ind,per)'));
CorrPermDistTest(C,D,rep*2,'ProteomicsPancvsNormalNeuroganin');

% Screening Panc vs Normal
pMGC=zeros(m,1);pD=zeros(m,1);pM=zeros(m,1);pP=zeros(m,1);pHHG=zeros(m,1);testMGC=zeros(m,1);testD=zeros(m,1);testM=zeros(m,1);testHHG=zeros(m,1);
for i=1:m
    i
    C=squareform(pdist(A(i,per)'));   
    [pMGC(i),pD(i),pM(i),pP(i), pHHG(i),testMGC(i),testD(i),testM(i),testHHG(i)]=CorrPermDistTest(C,D,rep);
end

testMGC(:,6)=pMGC;
testMGC(:,7)=pP;
testMGC(:,8)=pD;
testMGC(:,9)=pM;
testMGC(:,10)=pHHG;
testMGC(:,2)=testD;
testMGC(:,3)=testD;
testMGC(:,4)=testM;
testMGC(:,5)=testHHG;
save('ScreeningPancvsNormal');


clear
load proteomics.mat
% Panc vs All
per=(LabelIndAll~=2) & (LabelIndAll<5);
LabelIndAll(per)=1;
per=(LabelIndAll<5);
D=LabelIndAll(per);
D=squareform(pdist(D));
rep=5000;
m=318;

% C=squareform(pdist(A(:,per)'));
% CorrPermDistTest(C,D,rep*2,'ProteomicsOvarvsAll');
ind=181;
% C=A(:,per)';
% for i=1:318
% C(:,i)=C(:,i)/norm(C(:,i));
% end
C=squareform(pdist(A(ind,per)'));
CorrPermDistTest(C,D,rep*2,'ProteomicsPancvsAllNeuroganin');

% Screening Panc vs All
pMGC=zeros(m,1);pD=zeros(m,1);pM=zeros(m,1);pP=zeros(m,1);pHHG=zeros(m,1);testMGC=zeros(m,1);testD=zeros(m,1);testM=zeros(m,1);testHHG=zeros(m,1);
for i=1:m
    i
    C=squareform(pdist(A(i,per)'));   
    [pMGC(i),pD(i),pM(i),pP(i), pHHG(i),testMGC(i),testD(i),testM(i),testHHG(i)]=CorrPermDistTest(C,D,rep);
end

testMGC(:,6)=pMGC;
testMGC(:,7)=pP;
testMGC(:,8)=pD;
testMGC(:,9)=pM;
testMGC(:,10)=pHHG;
testMGC(:,2)=testD;
testMGC(:,3)=testD;
testMGC(:,4)=testM;
testMGC(:,5)=testHHG;
save('ScreeningPancvsAll');


thres=0.05;
min=1/rep;
seq=318:-1:1;
thres=thres./seq;
thres(thres<min)=min;
num=zeros(5,1);
for i=1:5
    tmp=testMGC(:,5+i);
    [tmp,ind]=sort(tmp,'ascend');
    num(i)=find(tmp-thres>0,1,'first')-1;
    ind=ind(1:num(i));
    num(i)=length(ind);
    
%     tmp=testMGC2(:,5+i);
%     [tmp,ind2]=sort(tmp,'ascend');
%     tmpNum=find(tmp-thres>0,1,'first')-1;
%     ind2=ind2(1:tmpNum);
%     ind2=intersect(ind,ind2);
%     num(i)=length(ind2);
end








load('BrainHippoShape.mat')
% [A,B,RX,RY]=MGCDistTransform(LMRS,LMLS);
% k=60;
% LMRS(RX>k)=0;
per1=(Label==1);
per2=(Label==2);
per3=(Label==3);
t1=LMRS(per1,per1);
t2=LMRS(per1,per2);
t3=LMRS(per1,per3);
t1=t1(t1>0);t2=t2(t2>0);t3=t3(t3>0);
t2=[t2;t3];
% a1=t1;a2=t2;a3=t3;
t1=LMRS(per2,per2);
t2=LMRS(per2,per1);
t3=LMRS(per2,per3);
t1=t1(t1>0);t2=t2(t2>0);t3=t3(t3>0);
t2=[t2;t3];
% a1=[a1;t1];a2=[a2;t2];a3=[a3;t3];
t1=LMRS(per3,per3);
t2=LMRS(per3,per2);
t3=LMRS(per3,per1);
t1=t1(t1>0);t2=t2(t2>0);t3=t3(t3>0);
t2=[t2;t3];
% a1=[a1;t1];a2=[a2;t2];a3=[a3;t3];

thres=500;
mean(t1>thres)
mean(t2>thres)
% mean(a3>thres)
% a1=[a1',a2']';
% hold on
% [f,xi]=ksdensity(t1);
% plot(xi,f,'r')
% [f,xi]=ksdensity(t2);
% plot(xi,f,'b')
% legend('class1 vs class 2','class2 vs class 3', 'class3 vs class1')
% [~,b,c]=kstest2(a1,a2)
% [~,b,c]=kstest2(a1,a3)
% [~,b,c]=kstest2(a2,a3)

pv=t1;
[f,xi]=ksdensity(pv);
hold on
h1=plot(xi,f,'b.-','LineWidth',2);
% pv=sort(pv,'ascend');
% ord=0.002*ones(length(pv),1);
% for i=2:length(pv)
%     if pv(i)-pv(i-1)<0.1
%         ord(i)=ord(i-1)+0.0001;
%     else
%         ord(i)=0.002;
%     end
% end
% plot(pv,ord,'b.','MarkerSize',2);

pv=t2;
[f,xi]=ksdensity(pv);
h2=plot(xi,f,'r.-','LineWidth',2);
% pv=sort(pv,'ascend');
% ord=0.001*ones(length(pv),1);
% for i=2:length(pv)
%     if pv(i)-pv(i-1)<0.1
%         ord(i)=ord(i-1)+0.0001;
%     else
%         ord(i)=0.001;
%     end
% end
% plot(pv,ord,'r.','MarkerSize',2);

title('Within-Class vs Between-Class Distances of Class 2')
set(gca,'YTick',[0,0.002,0.004],'YTickLabel',[0,0.002,0.004],'fontSize',12)
% set(gca,'YTick',[-0.004,-0.002,0,0.002,0.004],'YTickLabel',[0.004,0.002,0,0.002,0.004],'fontSize',12)
legend([h1,h2],'Within-Class Distances','Between-Class Distances')
% <<<<<<< HEAD
% xlim([0,0.15]);
% ylim([-1 max(f)+1]);



n1=10;
x=unifrnd(-1,1,n1,1);
y=x;
n2=2;
x2=unifrnd(3,4,n2,1);
y2=x2;
X=squareform(pdist([x;x2]));
Y=squareform(pdist([y;y2]));
[A,B,RX,RY]=MGCDistTransform(X,Y,'mgc');
C=A.*B;
A=A.*A;
B=B.*B;
CJ=sum(C,1)'./sqrt(sum(A,1)'.*sum(B,2));
CI=sum(C,2)./sqrt(sum(A,1)'.*sum(B,2));
hold on
plot(x,y,'k.');
plot(x2,y2,'r.');

n1=10;
x=unifrnd(-1,1,n1,1);
y=x;
X=squareform(pdist(x));
Y=squareform(pdist(y));
[A,B,RX,RY]=MGCDistTransform(X,Y,'mgc');
C=A.*B;
[corr,varX,varY]=MGCLocalCorr(X,Y,'mgc');
max(max(corr))-corr(end,end);
cov=corr.*sqrt(varX'*varY);
max(max(cov))-cov(end,end)
% A=A.*A;
% B=B.*B;
% CJ=sum(C,1)'./sqrt(sum(A,1)'.*sum(B,2));
% CI=sum(C,2)./sqrt(sum(A,1)'.*sum(B,2));






n=100;
x=unifrnd(0,1,100);
x=unifrnd(0,1,100,1);
y=x;
x2=unifrnd(0,1,100,1);
C=squareform(pdist(x+x2));
C=squareform(pdist([x';x2']'));
D=squareform(pdist(y));
dcov=sum(sum(C.*D))/n^2;
dvar1=(sum(sum(C.*C))/n^2)^0.5;
dvar2=(sum(sum(D.*D))/n^2)^0.5;

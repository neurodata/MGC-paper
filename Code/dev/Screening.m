load proteomics.mat

per=(LabelIndAll==1 | LabelIndAll==4);
D=LabelIndAll(per);
D=squareform(pdist(D));
D(D>0)=1;
rep=1000;
m=318;

pMGC=zeros(m,1);pD=zeros(m,1);pM=zeros(m,1);pP=zeros(m,1);pHHG=zeros(m,1);testMGC=zeros(m,1);testD=zeros(m,1);testM=zeros(m,1);testHHG=zeros(m,1);
for i=1:m
    i
    C=squareform(pdist(A(i,per)'));   
    [pMGC(i),pD(i),pM(i),pP(i), pHHG(i),testMGC(i),testD(i),testM(i),testHHG(i)]=CorrPermDistTest(C,D,rep);
end

testMGC(:,2)=pMGC;
testMGC(:,3)=testM;
testMGC(:,4)=pM;
testMGC(:,6)=pD;
testMGC(:,8)=pHHG;
testMGC(:,5)=testD;
testMGC(:,7)=testHHG;

C=squareform(pdist(A(:,per)'));
CorrPermDistTest(C,D,rep,'Proteomics1vs4');



load('BrainHippoShape.mat')
[A,B,RX,RY]=MGCDistTransform(LMRS,LMLS);
% k=50;
% LMRS(RX>k)=0;
per1=(Label==1);
per2=(Label==2);
per3=(Label==3);
a1=LMRS(per1,per1);
a2=LMRS(per2,per2);
a3=LMRS(per3,per3);
a1=a1(a1>0);a2=a2(a2>0);a3=a3(a3>0);
a1=[a1',a2']';
hold on
[f,xi]=ksdensity(a1);
plot(xi,f,'r')
[f,xi]=ksdensity(a2);
plot(xi,f,'b')
[f,xi]=ksdensity(a3);
plot(xi,f,'k')
legend('class1 vs class 2','class2 vs class 3', 'class3 vs class1')
[~,b,c]=kstest2(a1,a2)
[~,b,c]=kstest2(a1,a3)
[~,b,c]=kstest2(a2,a3)

pv=a3;
[f,xi]=ksdensity(pv);
hold on
plot(xi,f,'.-','LineWidth',3);
pv=sort(pv,'ascend');
ord=0.001*ones(length(pv),1);
for i=2:length(pv);
    if pv(i)-pv(i-1)<0.001
        ord(i)=ord(i-1)+0.001;
    end
end
plot(pv,ord,'.','MarkerSize',2);
% <<<<<<< HEAD
xlim([0,0.15]);
ylim([-1 max(f)+1]);

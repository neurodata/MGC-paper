%Locally linear dependent
n=7;
x=-1:2/(n-1):1;
y=x.^2;
x=unifrnd(1,2,1,n);
x=[-x,x];
y=abs(x);
%x=unifrnd(-1,1,1,n);
%y=x;
optionModi=1;
C=squareform(pdist(x'));
P=squareform(pdist(y'));
per=perms(1:n);
t1=LocalGraphCorr(C, P, optionModi);
pv=zeros(n,n);
for i=1:size(per,1);
    t2=LocalGraphCorr(C, P(per(i,:),per(i,:)), optionModi);
    pv=pv+(t1>t2);
end
pv=1-pv/size(per,1);
%%%ind trial
clear
load('BrainCP') 
n=42;
alpha=0.05; rep2=100;repp=100;
dim=1;type=20;
power1=0;power2=0;power3=0;power4=0;power5=0;power6=0;power7=0;power8=0;
for i=1:rep2;
    %x=random('norm',0,1,n,dim);
    %y=random('norm',0,1,n,dim);
    %y=x.^2;
    %[x, y]=CorrSampleGenerator(type,n,dim,1,0);
    y=mvnrnd(zeros(n,dim),eye(dim),n);
    distPInd=squareform(pdist(y));
    [p1, p2, p3, p4,n1,n2,n3]=CorrPermDistTest([distC distPInd],repp,repp,'BrainCxPInd');
%     distC=squareform(pdist(x));
%     distP=squareform(pdist(y));
%     [p1, p2, p3, p4,n1,n2,n3]=CorrPermDistTest([distC distP],rep2,repp);
    %if i==1
        neighbor1=n1;
        neighbor2=n2;
        neighbor3=n3;
    %end
    if mean(p1(neighbor1))<alpha
        power1=power1+1/rep2;
    end
    if mean(p2(neighbor2))<alpha
        power2=power2+1/rep2;
    end
    if mean(p3(neighbor3))<alpha
        power3=power3+1/rep2;
    end
    if p1(end,end)<alpha
        power4=power4+1/rep2;
    end
    if p2(end,end)<alpha
        power5=power5+1/rep2;
    end
    if p3(end,end)<alpha
        power6=power6+1/rep2;
    end
    if p4<alpha
        power7=power7+1/rep2;
    end
    if min(min(p2))<alpha
        power8=power8+1/rep2;
    end
end

%%Sims
%Ind
n=100; dim=1; lim=20; rep1=2000;rep2=10000;
CorrIndTest(1,n,dim,lim,rep1,rep2);
CorrIndTest(2,n,dim,lim,rep1,rep2);
CorrIndTest(3,n,dim,lim,rep1,rep2);
CorrIndTest(4,n,dim,lim,rep1,rep2);
CorrIndTest(5,n,dim,lim,rep1,rep2);
CorrIndTest(6,n,dim,lim,rep1,rep2);
CorrIndTest(7,n,dim,lim,rep1,rep2);
CorrIndTest(8,n,dim,lim,rep1,rep2);
CorrIndTest(9,n,dim,lim,rep1,rep2);
CorrIndTest(10,n,dim,lim,rep1,rep2);

n=100; dim=1; lim=20; rep1=2000;rep2=10000;
CorrIndTest(11,n,dim,lim,rep1,rep2);
CorrIndTest(12,n,dim,lim,rep1,rep2);
CorrIndTest(13,n,dim,lim,rep1,rep2);
CorrIndTest(14,n,dim,lim,rep1,rep2);
CorrIndTest(15,n,dim,lim,rep1,rep2);
CorrIndTest(16,n,dim,lim,rep1,rep2);
CorrIndTest(17,n,dim,lim,rep1,rep2);
CorrIndTest(18,n,dim,lim,rep1,rep2);
CorrIndTest(19,n,dim,lim,rep1,rep2);
CorrIndTest(20,n,dim,lim,rep1,rep2);



%IndDim
n=100; dim=1000; lim=20; rep1=2000;rep2=10000;
CorrIndTestDim(1,n,dim,lim,rep1,rep2);
CorrIndTestDim(2,n,dim,lim,rep1,rep2);
CorrIndTestDim(3,n,dim,lim,rep1,rep2);
CorrIndTestDim(4,n,dim,lim,rep1,rep2);
dim=10;lim=10;
CorrIndTestDim(5,n,dim,lim,rep1,rep2);
dim=40;lim=20;
CorrIndTestDim(6,n,dim,lim,rep1,rep2);
CorrIndTestDim(7,n,dim,lim,rep1,rep2);
CorrIndTestDim(8,n,dim,lim,rep1,rep2);
dim=40;lim=20;
CorrIndTestDim(9,n,dim,lim,rep1,rep2);
dim=100;
CorrIndTestDim(10,n,dim,lim,rep1,rep2);

n=100; dim=100; lim=20; rep1=2000;rep2=10000;
CorrIndTestDim(11,n,dim,lim,rep1,rep2);
dim=20;
CorrIndTestDim(12,n,dim,lim,rep1,rep2);
dim=10;lim=10;
CorrIndTestDim(13,n,dim,lim,rep1,rep2);
dim=20;lim=20;
CorrIndTestDim(14,n,dim,lim,rep1,rep2);
CorrIndTestDim(15,n,dim,lim,rep1,rep2);
CorrIndTestDim(16,n,dim,lim,rep1,rep2);
dim=20;
CorrIndTestDim(17,n,dim,lim,rep1,rep2);
dim=100;
CorrIndTestDim(18,n,dim,lim,rep1,rep2);
dim=100;
CorrIndTestDim(19,n,dim,lim,rep1,rep2);
CorrIndTestDim(20,n,dim,lim,rep1,rep2);


%%%Reals
clear
load('BrainBNU1.mat') 
C1=Corr1(:,:,1);
P1=Corr2(:,:,1);
P2=Corr2(:,:,2);rep2=50;
CorrPermDistTest([C1 P1],rep2,rep2,'BNU111');
CorrPermDistTest([C1 P2],rep2,rep2,'BNU112');

%%%%
clear
load('ccidiff-Tmat-org')
n=109; lim=1; rep1=1; rep2=2000;
CorrPermDistTest(C1,T1,rep2,rep2,'CT');

%%%use dcorr to find the optimal neighborhood size
clear
load('BrainCP') 
n=42; lim=1; rep1=2000; rep2=10000;
CorrPermDistTest(distC,distP,rep1,rep2,'BrainCxP');
% mean(p1(neighbor1))
% mean(p2(neighbor2))
%%%ind trial
d=1;
y=mvnrnd(zeros(n,d),eye(d),n);
distPInd=squareform(pdist(y));
CorrPermDistTest([distC distPInd],rep1,rep2,'BrainCxPInd');

%%%
clear
load('BrainHippoShape')
n=114;lim=1; rep1=2000;rep2=10000;
y=squareform(pdist(Label));
%y=(y>0)+1;
y=y+1;
for i=1:n
    y(i,i)=0;
end
yind=unifrnd(0,3,n,1);
yind=squareform(pdist(ceil(yind)));
%yind=(yind>0);
CorrPermDistTest(LMLS,y,rep1,rep2,'BrainLMLxY');
CorrPermDistTest(LMRS,y,rep1,rep2, 'BrainLMRxY');
CorrPermDistTest([LMLS LMRS],rep1,rep2,'BrainLMLxLMR');
%%%ind trial
CorrPermDistTest([LMLS yind], rep1,rep2,'BrainLMLxYInd');
CorrPermDistTest([LMRS yind], rep1,rep2,'BrainLMRxYInd');
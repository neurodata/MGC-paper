%Locally linear dependent
n=7;
x=-1:2/(n-1):1;
% x=unifrnd(-1,1,1,n);
% y=x.^2;
% x=unifrnd(1,2,1,n);
% x=[-x,x];
% y=abs(x);
%x=unifrnd(-1,1,1,n);
y=x;
optionModi=2;
C=squareform(pdist(x'));
P=squareform(pdist(y'));
per=perms(1:n);
t1=LocalCorr(C, P, optionModi);
pv=zeros(size(t1));
for i=1:size(per,1);
    t2=LocalCorr(C, P(per(i,:),per(i,:)), optionModi);
    pv=pv+(t1>t2);
end
pv=1-pv/size(per,1);
rat(pv)
%%%ind trial
clear
load('BrainCP')
n=30;
alpha=0.05; rep1=100;rep2=300;
dim=1;type=20;
power=zeros(7,1);
p=zeros(7,1);
%power=zeros(n,n,length(ratio));
for i=1:rep1;
    x=random('norm',0,1,n,dim);
    y=random('norm',0,1,n,dim);
    z=random('norm',0,1,n,dim)/n;
    %x=unifrnd(-1,1,n,1);
    %y=x.^2;
    %     %[x, y]=CorrSampleGenerator(type,n,dim,1,0);
    %     %y=mvnrnd(zeros(n,dim),eye(dim),n);
    distC=squareform(pdist(x));
    distPInd=squareform(pdist(y));
    [p(1),p(2),p(3),p(4),p(5),p(6),p(7)]=CorrPermDistTest3(distC,distPInd,5,rep2,'BrainCxPInd');
    for j=1:7
        if p(j)<alpha
            power(j)=power(j)+1/rep1;
            if j==1
                p(1)
                p(4)
            end
        end
    end
end

%%Sims
%outlier model
n=100; dim=1; lim=1; rep1=200;rep2=1000;
noise=0.3;
CorrIndTest(0,n,dim,lim,rep1,rep2,noise);
noise=0.5;
CorrIndTest(0,n,dim,lim,rep1,rep2,noise);
noise=0.7;
CorrIndTest(0,n,dim,lim,rep1,rep2,noise);
%Ind
n=100; dim=1; lim=20; rep1=200;rep2=1000;
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

n=100; dim=1; lim=20; rep1=200;rep2=1000;
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
n=100; dim=1000; lim=20; rep1=200;rep2=1000;
CorrIndTestDim(1,n,dim,lim,rep1,rep2);
CorrIndTestDim(2,n,dim,lim,rep1,rep2);
CorrIndTestDim(3,n,dim,lim,rep1,rep2);
dim=20;lim=20;
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

n=100; dim=20; lim=20; rep1=200;rep2=1000;
CorrIndTestDim(11,n,dim,lim,rep1,rep2);
CorrIndTestDim(12,n,dim,lim,rep1,rep2);
CorrIndTestDim(13,n,dim,lim,rep1,rep2);
dim=40; lim=20;
CorrIndTestDim(14,n,dim,lim,rep1,rep2);
CorrIndTestDim(15,n,dim,lim,rep1,rep2);
dim=10;lim=10;
CorrIndTestDim(16,n,dim,lim,rep1,rep2);
dim=10;lim=10;
CorrIndTestDim(17,n,dim,lim,rep1,rep2);
dim=1000;lim=20;
CorrIndTestDim(18,n,dim,lim,rep1,rep2);
dim=100;
CorrIndTestDim(19,n,dim,lim,rep1,rep2);
CorrIndTestDim(20,n,dim,lim,rep1,rep2);

%%
thres=0.8;dim=1;rep1=100;rep2=200;noise=1;type=1:10;
[p1]=CorrSimPermScale1(type,dim,thres,rep1,rep2,noise);
thres=0.8;dim=1;rep1=100;rep2=200;noise=1;type=11:20;
[p1]=CorrSimPermScale1(type,dim,thres,rep1,rep2,noise);
thres=0.5;dim=2;rep1=100;rep2=200;noise=0;type=1:10;
[p1]=CorrSimPermScale1(type,dim,thres,rep1,rep2,noise);
thres=0.5;dim=2;rep1=100;rep2=200;noise=0;type=11:20;
[p1]=CorrSimPermScale1(type,dim,thres,rep1,rep2,noise);

%%%use dcorr to find the optimal neighborhood size
clear
load('BrainCP')
n=42; rep=10000;
CorrPermDistTest(distC,distP,rep,'BrainCxP');
% mean(p1(neighbor1))
% mean(p2(neighbor2))

%%%
clear
load('semipar')
n=109;rep=10000;
distCCI=squareform(pdist(cci));
CorrPermDistTest(distMigrain(ind,ind),distCCI(ind,ind),rep,'MigrainxCCI');
CorrPermDistTest(distM2g(ind,ind),distCCI(ind,ind),rep,'M2gxCCI');
CorrPermDistTest(distM2g(ind,ind),distMigrain(ind,ind),rep,'M2gxMigrain');

%%%
clear
load('BrainHippoShape')
n=114;rep=10000;alpha=0.05;
y=squareform(pdist(Label));
% y=(y>0)+1;
y=y+1;
for i=1:n
    y(i,i)=0;
end
% y(y>0)=1;
%estimate optimal scale separately
option=[1,2,3,4];
CorrPermDistTest(LMLS,y,rep,'BrainLMLxY',option);
CorrPermDistTest(LMRS,y,rep, 'BrainLMRxY',option);
CorrPermDistTest(LMLS,LMRS,rep,'BrainLMLxLMR');
%%%ind trial
rep1=200;rep2=200;powerL=zeros(7,1);powerR=zeros(7,1);
for i=1:rep2;
    yind=unifrnd(0,3,n,1);
    yind=squareform(pdist(ceil(yind)));
    [p(1),p(2),p(3),p(4),p(5),p(6),p(7)]=CorrPermDistTest(LMRS,yind,rep1,rep2,'BrainLMRxYIndJ');
    for j=1:7
        if p(j)<alpha
            powerR(j)=powerR(j)+1/rep2;
        end
    end
    
    [p(1),p(2),p(3),p(4),p(5),p(6),p(7)]=CorrPermDistTest(LMLS,yind,rep1,rep2,'BrainLMLxYIndJ');
    for j=1:7
        if p(j)<alpha
            powerL(j)=powerL(j)+1/rep2;
        end
    end
end


%%%knn trial
clear
load('BrainHippoShape')
n=114;lim=1; rep1=2000;rep2=10000;alpha=0.05;
error=CorrKNN(LMLS,Label);



load('Wiki_Data')
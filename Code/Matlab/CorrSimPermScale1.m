function [powerP]=CorrSimPermScale1(n,dim,type,rep1,rep2,noise,alpha)
% Author: Cencheng Shen
% n=50;dim=1;rep1=50;rep2=300;noise=1;type=1:20;
% [p1]=CorrSimPermScale1(n,dim,type,rep1,rep2,noise);
% n=100;dim=10;rep1=50;rep2=300;noise=0;type=1:20;
% [p1]=CorrSimPermScale1(n,dim,type,rep1,rep2,noise);
% Used to compute the permutation test power for the 20 type of simulations
if nargin<1
    n=50;
end
if nargin<2
    dim=1;
end
if nargin<3
    type=1:20;
end
if nargin<4
    rep1=50;
end
if nargin<5
    rep2=300; 
end
if nargin<6
    noise=1;
end
if nargin<7
    alpha=0.05; % Default type 1 error level
end
if dim>1
    noise=0;
    dimInd=dim;
end
repp=4;
pre1='../../Data/'; 
%pre2='../../Figures/Fig'; % The folder to save figures
powerP=zeros(6,20);
option=[1,2,0,0];
for tt=type
    if dim==1
        filename=strcat(pre1,'CorrIndTestType',num2str(tt),'N',num2str(100),'Dim',num2str(dim));
        load(filename,'neighborhoods','numRange','power1All','power2All','power3All','power4','power5','power6','power7');
        ind=find(numRange==n);
        neighbor=neighborhoods(:,ind);
        power1All=power1All(1:numRange(ind),1:numRange(ind),ind);
        power2All=power2All(1:numRange(ind),1:numRange(ind),ind);
        power3All=power3All(1:numRange(ind),1:numRange(ind),ind);
    else
         filename=strcat(pre1,'CorrIndTestDimType',num2str(tt),'N',num2str(100),'Dim');
         load(filename,'neighborhoods','dimRange','power1All','power2All','power3All','power4','power5','power6','power7');
         if dimInd>length(dimRange)
             dimInd=length(dimRange);
         end
         neighbor=neighborhoods(:,dimInd);
         dim=dimRange(dimInd);
         power1All=power1All(:,:,dimInd);
         power2All=power2All(:,:,dimInd);
         power3All=power3All(:,:,dimInd);
         ind=dimInd;
    end
    p=zeros(6,1);
    tt
    %CorrIndTest(tt,n,dim,1,rep2,0);
    for r1=1:rep1
        [x, y]=CorrSampleGenerator(tt,n,dim,1, noise);
        C=squareform(pdist(x));
        D=squareform(pdist(y));
       [pp1,pp2,~,p1,p2]=CorrPermDistTest(C,D,repp,rep2,'PermInd',option);
%          [pp1,~,~,p2]=CorrPermDistTest(C,D,0,rep2,'PermInd',ratio,option);
        p(1)=p(1)+(pp1<alpha)/rep1;
        p(2)=p(2)+(p1<alpha)/rep1;
        p(5)=p(5)+(pp2<alpha)/rep1;
        p(6)=p(6)+(p2<alpha)/rep1;
    end
    powerP(1,tt)=p(1);
    powerP(2,tt)=p(2);
    powerP(5,tt)=p(5);
    powerP(6,tt)=p(6);
    if n<=100
        powerP(3,tt)=max(max(power1All));
        powerP(4,tt)=power1All(end,end);
    end
end
filename=strcat(pre1,'CorrSimPermScale',num2str(type(1)),'-',num2str(type(end)),'N',num2str(n),'Dim',num2str(dim));
save(filename,'powerP','n','rep1','rep2','dim','noise','alpha');

figure
x=1:20;
p1=powerP;
plot(x,p1(1,:),'bo-',x,p1(2,:),'rx--',x,p1(3,:),'ko-',x,p1(4,:),'cx--')
legend('Estimated MGC{dcorr}','Estimated dCorr', 'True MGC{dcorr}','True dCorr');
ylim([0,1]);
title('Testing Power Comparison for dimension 1 simulation at n=50');

%  k=n;l=n;
% load('tmpInd.mat')
% figure
% ksdensity(reshape(dCor1A(k,l,:),1,length(dCor1A(k,l,:))))
% hold on
% ksdensity(reshape(dCor1N(k,l,:),1,length(dCor1A(k,l,:))))
% xlim([-0.1 1])
% n1=maxNeighbors(power1,dCor1N,dCor1A)
% load('tmpPerm.mat')
% figure
% hold on
% ksdensity(reshape(dCor1A(k,l,:),1,length(dCor1A(k,l,:))))
% ksdensity(reshape(dCor1N(k,l,:),1,length(dCor1A(k,l,:))))
% xlim([-0.1 1])
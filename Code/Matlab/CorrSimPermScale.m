function [powerP,powerI]=CorrSimPermScale(n,dim,type,rep1,rep2,noise,alpha)
% Author: Cencheng Shen
% n=50;dim=1;rep1=100;rep2=300;noise=1;type=1:5;
% [p1,p2]=CorrSimPermScale(n,dim,type,rep1,rep2,noise);
% n=50;dim=1;rep1=100;rep2=300;noise=1;type=6:10;
% [p1,p2]=CorrSimPermScale(n,dim,type,rep1,rep2,noise);
% n=50;dim=1;rep1=100;rep2=300;noise=0;type=11:15;
% [p1,p2]=CorrSimPermScale(n,dim,type,rep1,rep2,noise);
% n=50;dim=1;rep1=100;rep2=300;noise=0;type=16:20;
% [p1,p2]=CorrSimPermScale(n,dim,type,rep1,rep2,noise);
% n=30;dim=10;rep1=100;rep2=300;noise=0;type=1:19;
% [p1]=CorrSimPermScale(n,dim,type,rep1,rep2,noise);
% Used to compute the permutation test power for the 20 type of simulations
if nargin<1
    n=100;
end
if nargin<2
    dim=1;
end
if nargin<3
    type=1:20;
end
if nargin<4
    rep1=100;
end
if nargin<5
    rep2=100; 
end
if nargin<6
    noise=1;
end
if nargin<7
    alpha=0.05; % Default type 1 error level
end

total=20;
if dim>1
    noise=0;
end

pre1='../../Data/'; 
%pre2='../../Figures/Fig'; % The folder to save figures
powerP=zeros(7,total);
powerI=zeros(7,total);
option=[1,1,1,1,1,1,1];
for tt=type
    filename=strcat(pre1,'CorrIndTestType',num2str(tt),'N',num2str(100),'Dim',num2str(dim));
    load(filename,'neighborhoods','numRange');
    neighbor=neighborhoods(:,find(numRange==n));
    p=zeros(7,1);
    tt
    %CorrIndTest(tt,n,dim,1,rep2,0);
    for r1=1:rep1
        r1
        [x, y]=CorrSampleGenerator(tt,n,dim,1, noise);
        [p(1),p(2),p(3),p(4),p(5),p(6),p(7)]=CorrPermDistTest(squareform(pdist(x)),squareform(pdist(y)),rep2,rep2,'PermInd',alpha,option);
        for j=1:7
            if p(j)<alpha
                powerP(j,tt)=powerP(j,tt)+1/rep1;
            end
        end
        
        [p(1),p(2),p(3),p(4),p(5),p(6),p(7)]=CorrPermDistTest(squareform(pdist(x)),squareform(pdist(y)),rep2,rep2,'PermInd',alpha,option,neighbor);
        for j=1:7
            if p(j)<alpha
                powerI(j,tt)=powerI(j,tt)+1/rep1;
            end
        end
    end
end

filename=strcat(pre1,'CorrSimPermScale',num2str(type(1)),'-',num2str(type(end)),'N',num2str(n),'Dim',num2str(dim));
save(filename,'powerP','powerI','n','rep1','rep2','dim','noise','alpha');

%  k=n;l=n;
% % load('tmpInd.mat')
% % figure
% % ksdensity(reshape(dCor1A(k,l,:),1,length(dCor1A(k,l,:))))
% % hold on
% % ksdensity(reshape(dCor1N(k,l,:),1,length(dCor1A(k,l,:))))
% % xlim([-0.1 0.5])
% % n1=maxNeighbors(power1,dCor1N,dCor1A)
% load('tmpPerm.mat')
% figure
% hold on
% ksdensity(reshape(dCor1A(k,l,:),1,length(dCor1A(k,l,:))))
% ksdensity(reshape(dCor1N(k,l,:),1,length(dCor1A(k,l,:))))
%xlim([-0.1 0.5])
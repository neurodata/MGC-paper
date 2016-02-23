function [power]=CorrSimPermScale(n,dim,type,rep1,rep2,noise,alpha)
% Author: Cencheng Shen
% n=100;dim=1;rep1=100;rep2=100;noise=1;type=1:1;
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
power=zeros(7,total);
for tt=type
    p=zeros(7,1);
    tt
    for r1=1:rep1
        [x, y]=CorrSampleGenerator(tt,n,dim,1, noise);
        [p(1),p(2),p(3),p(4),p(5),p(6),p(7)]=CorrPermDistTest(squareform(pdist(x)),squareform(pdist(y)),rep2,rep2,'PermInd',alpha);
        for j=1:7
            if p(j)<alpha
                power(j,tt)=power(j,tt)+1/rep1;
            end
        end
    end
end

filename=strcat(pre1,'CorrSimPermScale',num2str(length(type)),'N',num2str(n),'Dim',num2str(dim));
save(filename,'power','n','rep1','rep2','dim','noise','alpha');
function [powerP]=CorrSimPermTest(type,dim,thres,rep1,rep2,noise,alpha)
% Author: Cencheng Shen
% n=60;dim=1;rep1=100;rep2=200;noise=1;type=1:20;
% [p1]=CorrSimPermScale(n,dim,type,rep1,rep2,noise);
% n=100;dim=10;rep1=100;rep2=200;noise=0;type=1:20;
% [p1]=CorrSimPermScale(n,dim,type,rep1,rep2,noise);
% Used to compute the permutation test power for the 20 type of simulations
if nargin<1
    type=1:20;
end
if nargin<2
    dim=1;
end
if nargin<3
    thres=0.8;
end
if nargin<4
    rep1=1000; %1000
end
if nargin<5
    rep2=1000; %1000
end
if nargin<6
    noise=1;
end
if nargin<7
    alpha=0.05; % Default type 1 error level
end
if dim>1
    noise=0;
    n=100;
end
%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-3));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

pre1=strcat(rootDir,'Data/Results/');% The folder to save figures
powerP=zeros(7,20);
option=[0,2,0,0];
for tt=type
    neighbor=[];
    if dim==1
        filename=strcat(pre1,'CorrIndTestType',num2str(tt),'N',num2str(100),'Dim',num2str(dim));
        load(filename,'neighborhoods','numRange','power1All','power2All','power3All','power1','power2','power3','power7');
        ind=[find(max(power2,[],1)>=thres,1) length(numRange)];
        ind=min(ind);
        n=numRange(ind);
        neighbor=neighborhoods(:,ind);
        power1All=power1All(1:numRange(ind),1:numRange(ind),ind);
        power2All=power2All(1:numRange(ind),1:numRange(ind),ind);
        power3All=power3All(1:numRange(ind),1:numRange(ind),ind);
    else
        filename=strcat(pre1,'CorrIndTestDimType',num2str(tt),'N',num2str(100),'Dim');
        load(filename,'neighborhoods','dimRange','power1All','power2All','power3All','power1','power2','power3','power7','lim');
        ind=[find(power2>=thres,1,'last') 1];
        ind=max(ind);
        dim=dimRange(ind);
        neighbor=neighborhoods(:,ind);
        power1All=power1All(:,:,ind);
        power2All=power2All(:,:,ind);
        power3All=power3All(:,:,ind);
    end
    p=zeros(7,1);
    tt
    %CorrIndTest(tt,n,dim,1,rep2,0);
    for r1=1:rep1
        [x, y]=CorrSampleGenerator(tt,n,dim,1, noise);
        C=squareform(pdist(x));
        D=squareform(pdist(y));
        [pp1,pp2,~,p1,p2,~,p4,p1All,p2All]=CorrPermDistTest(C,D,rep2,'PermInd',option);
        p(1)=p(1)+(pp1<alpha)/rep1;
        p(2)=p(2)+(p1<alpha)/rep1;
        p(4)=p(4)+(pp2<alpha)/rep1;
        p(5)=p(5)+(p2<alpha)/rep1;
        if isempty(neighbor)==false
            if option(1)==1
                [k,l]=ind2sub(size(power1All),neighbor(1));
                k=min(k,size(p1All,1));l=min(l,size(p1All,2));
                p(3)=p(3)+(p1All(k,l)<alpha)/rep1;
            end
            if option(2)==2
                [k,l]=ind2sub(size(power2All),neighbor(2));
                k=min(k,size(p2All,1));l=min(l,size(p2All,2));
                p(6)=p(6)+(p2All(k,l)<alpha)/rep1;
            end
            if option(4)==4
                p(7)=p(7)+(p4<alpha)/rep1;
            end
        end
%         p
    end
    powerP(1,tt)=p(1);
    powerP(2,tt)=p(2);
    powerP(4,tt)=p(4);
    powerP(5,tt)=p(5);
    powerP(7,tt)=p(7);
    if isempty(neighbor)==false
        if option(1)==1
            powerP(3,tt)=p(3);
        end
        if option(2)==2
            powerP(6,tt)=p(6);
        end
    end
end
powerP(:,type)

filename=strcat(pre1,'CorrSimPermScale',num2str(type(1)),'-',num2str(type(end)),'Dim',num2str(dim));
save(filename,'powerP','rep1','rep2','dim','noise','alpha','thres');

figure
x=1:20;
a=2;
if a==1
    ind=[1,2,3,7];
else
    ind=[4,5,6,7];
end
p1=powerP(ind,:);
plot(x,p1(1,:),'bo-',x,p1(2,:),'rx--',x,p1(3,:),'ko-',x,p1(4,:),'c.-')
legend('Estimated MGC','Global Mcorr', 'True MGC','HHG');
xlabel('Function Type');
ylabel('Testing Power');
ylim([0,1]);
title('Testing Power Comparison for dimension 1 simulation at n=60');
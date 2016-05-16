function [ind]=MGCScaleVerify2(V)
% An auxiliary function to verify and estimate the MGC optimal scale
% cc=mean(mean(pAll(2:end,2:end)<0.05));
% if cc>0.05
%     ind=find(pAll==min(min(pAll)),1,'last');
% else
%     ind=1;
% end

VN=V(2:end,2:end);
VNN=VN(VN<0);

if (isempty(VNN))
    k=2:size(V,1);
    l=2:size(V,2);
else
    stdN=norm(VNN,'fro')/sqrt(length(VNN));
    VN=VN/stdN;
    k=VerifyRow(VN)+1;
    l=VerifyRow(VN')+1;
% %     [k,l]=size(V);
figure
imagesc(VN);
caxis([0 4])
colorbar();
end
ind=find(V==max(max(V(k,l))),1,'last');

function k=VerifyRow(tmp)
thres1=4;
thres2=0.5;

% t1=mean(tmp>thres1,2);
% if max(t1)<thres2  || sum(t1>thres2)<2
%     k=size(tmp,1);
% else
%     k=find(t1>thres2);
% end

[m,n]=size(tmp);
t=[];
for i=1:m;
    cE=find(tmp(i,:)<thres1);
    cS=[0 cE];
    cE=[cE m+1];
    c=max(cE-cS)-1;
    if c/n>thres2
        t=[t,i];
    end
end

k=[];
if (length(t)>1)
    if t(1)+1==t(2)
        k=[k t(1)];
    end
    for i=2:length(t)-1
        if (t(i)-1==t(i-1) && (t(i)+1==t(i+1)));
            k=[k t(i)];
        end
    end
    if t(length(t)-1)+1==t(length(t))
        k=[k t(end)];
    end
end

if (length(k)/n<0.1)
    k=m;
end

function [k,t]=VerifyRow3(pRow)
ll=length(pRow);
thres=0.6;
thres2=0.2;

% cE=find(pRow<thres);
% if size(cE,1)>1
%     cE=cE';
% end
% cS=[0 cE];
% cE=[cE ll+1];
% c=max(cE-cS)-1;
% t=c/ll;

t=length(find(pRow>thres))/ll;
% figure
% hist(pRow)
if t<thres2
    k=ll;
else
    k=find(pRow>thres);
end
% if k(end)==ll
%     k=ll;
% end

function [k,l,c,t]=VerifyRow2(pAll)
thres=0.85;
alpha=0.05;
t=0;
nX=size(pAll,1);
pRow=zeros(nX,1);
for i=1:nX
    pRow(i)=mean(pAll(i,2:end)<alpha);
end
% cE=find(pRow<thres)';
% cS=[0 cE];
% cE=[cE nX+1];
% c=max(cE-cS)-1;
% c=c/nX;
c=mean(pRow>thres);
if (pRow(end)>thres ||  sum(pRow(end)<pRow)==0)
    t=1;
end
k=find(pRow==max(pRow));
ind=find(pAll(k,:)==min(min(pAll(k,:))),1,'last');
[k1,l]=ind2sub(size(pAll(k,:)),ind);
k=k(k1);
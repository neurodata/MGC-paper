function [ind]=MGCScaleVerify(pAll)
% An auxiliary function to verify and estimate the MGC optimal scale
[k,l,c,t]=VerifyRow(pAll);
[lt,kt,ct,tt]=VerifyRow(pAll');
if c<ct
    l=lt;
    k=kt;
end
if c+ct<0.05
    k=1;
    l=1;
    if t==1 || tt==1
        [k,l]=size(pAll);
    end
end
% if c+ct>0
%     c+ct
% end
ind=sub2ind(size(pAll),k,l);

function [k,l,c,t]=VerifyRow(pAll)
thres=0.8;
alpha=0.05;
t=0;
nX=size(pAll,1);
pRow=zeros(nX,1);
for i=1:nX
    pRow(i)=mean(pAll(i,2:end)<alpha);
end
cE=find(pRow<thres)';
cS=[0 cE];
cE=[cE nX+1];
c=max(cE-cS)-1;
c=c/nX;
% c=mean(pRow>thres);
if (pRow(end)>thres ||  sum(pRow(end)<pRow)==0)
    t=1;
end
k=find(pRow==max(pRow));
ind=find(pAll(k,:)==min(min(pAll(k,:))),1,'last');
[k1,l]=ind2sub(size(pAll(k,:)),ind);
k=k(k1);
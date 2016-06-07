function [ind,pAll]=MGCScaleVerify2(V)
% An auxiliary function to verify and estimate the MGC optimal scale
VN=V(2:end,2:end);
stdN=(VN<=0);
tmp=sum(sum(stdN));
if tmp>0
    stdN=norm(VN(stdN))/sqrt(tmp);
else
    stdN=0.00001;
end

tmp=(VN>6*stdN);
% figure
% imagesc(tmp)
% mean(mean(tmp))
if mean(mean(tmp))>0.1%max(sum(tmp,1))> size(VN,1)/2 || max(sum(tmp,2))>size(VN,2)/2
    ind=find(V==max(max(VN)),1,'last');
else
    ind=size(V,1)*size(V,2);
end
pAll=ones(size(V));

[~,tmp]=sort(VN(VN<2),'descend');
ll=length(tmp);
ll2=size(VN);
for i=1:ll
    [k,l]=ind2sub(ll2,tmp(i));
    pAll(k+1,l+1)=(i-1)/(ll-1);
end


function k=Verify(VN)
thres=0.05;
[m,~]=size(VN);

k=m;

rowTmp=median(VN,2)';
indK=find(rowTmp==min(rowTmp),1,'last');
if rowTmp(indK)<=thres
    rowTmp=min(VN,[],2)';
    rowTmp=(rowTmp<thres);
    tmp=indK;
    for i=indK-1:-1:1
        if rowTmp(i)==1
            tmp=[i tmp];
        else
            break;
        end
    end
    for i=indK+1:m
        if rowTmp(i)==1
            tmp=[tmp i];
        else
            break;
        end
    end
    VN=VN(tmp,:);
    VN=VN(VN>-2);
    if median(VN)<=thres/m*length(tmp)
        k=tmp;
    end
end
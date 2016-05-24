function [ind]=MGCScaleVerify(V)
% An auxiliary function to verify and estimate the MGC optimal scale

VN=V(2:end,2:end);
k=Verify(VN)+1;
l=Verify(VN')+1;
ind=find(V==min(min(V(k,l))),1,'last');

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
function [ind]=MGCScaleVerify(V)
% An auxiliary function to verify and estimate the MGC optimal scale

VN=V(2:end,2:end);
k=Verify(VN)+1;
l=Verify(VN')+1;
ind=find(V==min(min(V(k,l))),1,'last');

function k=Verify(VN)
thres=0.05;
[m,~]=size(VN);
lim=2;

rowTmp=median(VN,2);
[~,indK]=sort(rowTmp,'ascend');

k=m;
for l=1:lim
    switch l
        case 1
            tmpThres=thres;
        case 2
            tmpThres=thres/2;
        case 3 
            tmpThres=thres/5;
    end
    len=indK(1:ceil(m/lim*(lim-l+1)));
    if median(rowTmp(len))<=tmpThres
        k=len;
        break;
    end
end
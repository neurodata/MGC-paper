function [p,indAll]=MGCScaleVerify(P)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations
thres=0.05;
thres2=0.005;

[m,n]=size(P);
tt=ceil(thres*[m,n]);
tt=[2,2];

p=P(end);
indAll=m*n;
% sum(sum(P<p))/(m-1)/(n-1)
if sum(sum(P<p))/(m-1)/(n-1)>thres
%     PN=floor(P*100)/100;
    
    C1=MonotoneRegion(P,thres2);
    C3=MonotoneRegion2(P,-thres2);
    if sum(sum(C3))>sum(sum(C1))
        C1=C3;
    end
    C2=MonotoneRegion(P',thres2)';
    C4=MonotoneRegion2(P',-thres2)';
    if sum(sum(C4))>sum(sum(C2))
        C2=C4;
    end
    if sum(sum(C2))>sum(sum(C1))
        C1=C2;
    end
    [~, ~, ~, C1] = FindLargestRectangles((C1==1), [0 0 1],tt);
%     figure
%     imagesc(C1)
    tmp=mean(mean(C1));
    if tmp >0
        p=prctile(P(C1), min(ceil(thres/tmp*100),100));
        [~, ~, ~, indAll]=FindLargestRectangles((C1==1) & (P<=p), [0 0 1]);
    end
end

function pdiff=MonotoneRegion(P,thres)
if nargin<2
    thres=0;
end
[m,n]=size(P);
pdiff=zeros(m,n);
for i=2:m
    tt=P(i,:);
    ttd=diff(tt);
    pdiff(i,2:end)=(ttd<=thres);
end

function pdiff=MonotoneRegion2(P,thres)
if nargin<2
    thres=0;
end
[m,n]=size(P);
pdiff=zeros(m,n);
for i=2:m
    tt=P(i,:);
    ttd=diff(tt);
    pdiff(i,2:end)=(ttd>=thres);
end
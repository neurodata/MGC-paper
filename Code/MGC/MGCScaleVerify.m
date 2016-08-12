function [p,indAll]=MGCScaleVerify(P)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations
thres=0.025;

[m,n]=size(P);
tt=ceil(thres*[m,n]);

p=min(max(max(P(2:end,2:end))),0.5);
indAll=1;

pc=mean(mean(P<=P(end)));
if pc<thres
    p=P(end);
    indAll=m*n;
else
    warning ('off','all');
    PN=floor(P*100)/100;
    
    C1=MonotoneRegion(PN,0);
    C2=MonotoneRegion(PN',0);
    if sum(sum(C2))>sum(sum(C1))
        C1=C2';
    end
    [~, ~, ~, C1] = FindLargestRectangles((C1==1), [0 0 1],tt);
    
    tmp=mean(mean(C1));
    if tmp >thres
        p=prctile(P(C1), ceil(thres/tmp*100));
        indAll=(C1==1) & (P<p);
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
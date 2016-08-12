function [p,indAll,t]=MGCScaleVerify(P)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations
thres=0.025;

[m,n]=size(P);
tt=ceil(thres*[m,n]);

p=min(max(max(P(2:end,2:end))),0.5);
% p=P(end);
indAll=1;

pc=mean(mean(P<=P(end)));
if pc<thres
    p=P(end);
    indAll=m*n;
else
%     prc=delta:delta:100-delta;
%     prc=[prc';pc;mean(mean(P<alpha))*100];
%     prc=sort(prc,'ascend');
%     prc=[mean(mean(P<alpha))*100];
    warning ('off','all');
    PN=floor(P*100)/100;
    
    C1=MonotoneRegion(PN,0);
    %[~, ~, ~, C1] = FindLargestRectangles((C1==1), [0 0 1],tt);
    C2=MonotoneRegion(PN',0);
    [~, ~, ~, C2] = FindLargestRectangles((C2==1), [0 0 1],tt);
    if sum(sum(C2))>sum(sum(C1))
        C1=C2';
    end
    
    tmp=mean(mean(C1));
    if tmp >thres
        p=prctile(P(C1), ceil(thres/tmp*100));
        indAll=(C1==1) & (P<p);
    end
%     for j=1:length(prc)
%         alpha=prctile(P(P<1), prc(j));
%         CC=(C1==1) & (P<=alpha);
% %         CC=(P<=alpha);
%          [~, ~, ~, CC] = FindLargestRectangles(CC, [0 0 1],ceil(thres*[m,n]));
% %         CC=BoundedRegion(CC,thres);
% %         CC=BoundedRegion(CC',thres)';
% %             try
% %                 CC=bwareafilt(CC,1);
% %             catch
% %                 CC=zeros(m,n);
% %             end
% %             figure
% %             imagesc(CC)
%         mm=sum(sum(CC))/(m-1)/(n-1);
%         if mm>thres
%             indAll=(CC==1);
%             p=alpha;
%             break;
%         end
%     end
end

% function CC2=BoundedRegion(CC,thres)
% [m,n]=size(CC);
% nn=ceil(thres*m);
% nn2=ceil(thres*n);
% CC2=(zeros(m,n)==1);
% i=2;
% while i<=m
%     is=max(2,i-nn);
%     ie=min(m,i+nn);
%     tmp=find(sum(CC(is:ie,:),1)<ie-is+1);
%     tmp=[1;tmp';n+1];
%     tmp2=diff(tmp);
%     ind=find(tmp2>=max(nn2*2));%,max(tmp2)));
%     i=i+1;
%     if ~isempty(ind)
%         for t=1:length(ind)
%             j=ind(t);
%             CC2(is:ie,tmp(j)+1:tmp(j+1)-1)=1;
%         end
%         i=ie+1;
%     end
% end

% 
function pdiff=MonotoneRegion(P,thres)
if nargin<2
    thres=0;
end
% thres=0.002;
[m,n]=size(P);
pdiff=zeros(m,n);
for i=2:m
    tt=P(i,:);
    ttd=diff(tt);
    pdiff(i,2:end)=(ttd<=thres); 
end
% figure
% imagesc(pdiff)
function [p,indAll]=MGCScaleVerify(V)
% An auxiliary function to verify and estimate the MGC optimal scale based
% on the p-values of all local correlations.
% VN=V(2:end,2:end);
% k=Verify(VN)+1; % verify the rows
% l=Verify(VN')+1; % verify the columns
% ind=find(V==min(min(V(k,l))),1,'last');

VN=V(2:end,2:end);
p1=Verify(VN);
p2=Verify(VN);
p=min(p1,p2);
%p=min(p,V(end));
indAll=find(V<=p,1,'last');
if isempty(indAll);
    indAll=size(V,1)*size(V,2);
end

function p=Verify(VN)
[m,~]=size(VN);
%k=m; % take the largest scale by default

rowTmp=median(VN,2);
rowTmp=sort(rowTmp,'ascend');
pp=ones(m,1);
tmp=0;
for i=1:m;
%     if mod(i,2) == 1
%         pp(i)=rowTmp(ceil(i/2));
%     else
%         pp(i)=(rowTmp(i/2)+rowTmp(i/2+1))/2;
%     end
%     pp(i)=pp(i)/i*m;
%     
    tmp=tmp+rowTmp(i);
    pp(i)=tmp/i/i*m;
end
p=min(pp);

% 
% function k=Verify(VN)
% thres=0.05; % type 1 error level
% [m,~]=size(VN);
% k=m; % take the largest scale by default
% 
% rowTmp=median(VN,2)';
% indK=find(rowTmp==min(rowTmp),1,'last');
% if rowTmp(indK)<=thres
%     rowTmp=min(VN,[],2)';
%     rowTmp=(rowTmp<thres);
%     tmp=indK;
%     
%     % from the row with minimal median p-value, include adjacency rows
%     % whose median p-value are significant
%     for i=indK-1:-1:1
%         if rowTmp(i)==1
%             tmp=[i tmp];
%         else
%             break;
%         end
%     end
%     for i=indK+1:m
%         if rowTmp(i)==1
%             tmp=[tmp i];
%         else
%             break;
%         end
%     end
%     
%     % check if the included rows satisfy a Bonferroni-type bound
%     VN=VN(tmp,:);
%     VN=VN(VN>-2);
%     if median(VN)<=thres/m*length(tmp)
%         k=tmp;
%     end
% end
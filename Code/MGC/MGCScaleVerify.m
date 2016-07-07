function [p,indAll]=MGCScaleVerify(V)
% V=enlarge(V);
% p=min(min(V(2:end,2:end)));
% indAll=find(V<=p);
% 
% function V=enlarge(V)
% [m,n]=size(V);
% V(1:m,1)=1;
% V(1,1:n)=1;
% for i=2:m
%     for j=2:n
%         V(i,j)=V(i,j)*max((m-1)/(i-1),(n-1)/(j-1))^2;
%     end
% end


% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations.
VN=V(2:end,2:end);
p1=Verify(VN); % verify the rows
p2=Verify(VN'); % verify the columns
p=min(p1,p2);
p=min(p,V(end));
indAll=find(V<=p);
if isempty(indAll)
    indAll=size(V,1)*size(V,2);
end
% p=V(indAll(end));

function p1=Verify(VN)
[m,~]=size(VN);

rowTmp=median(VN,2)';
for j=1:m
    rowTmp(j)=prctile(VN(j,:),55);
end
%rowTmp=min(VN,[],2)';
indK=find(rowTmp==min(rowTmp),1,'last');

% from the row with minimal median p-value, include adjacency rows
% whose median p-value are significant
tmpP=min(rowTmp)*m;
% ss=indK;
tt=indK;
for i=indK-1:-1:1
    tt=[i tt];
    tmp=mean(rowTmp(tt));
    if tmp*m/length(tt)<tmpP
        tmpP=tmp*m/length(tt);
%         ss=i;
    else
        break;
    end
end

% tmpP=min(rowTmp)*m;
% se=indK;
% tt=indK;
for i=indK+1:m
    tt=[tt i];
    tmp=mean(rowTmp(tt));
    if tmp*m/length(tt)<tmpP
        tmpP=tmp*m/length(tt);
        se=i;
    else
        break;
    end
end
% 
% indK=ss:1:se;
indK=tt;
p1=mean(rowTmp(indK))*m/length(indK);
% check if the included rows satisfy a Bonferroni-type bound
if p1>VN(end)
    p1=rowTmp(end);
end

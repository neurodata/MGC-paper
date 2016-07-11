function [p,indAll]=MGCScaleVerify(P,V)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations.
% m=mean(mean(VN));
% tmp=VN(VN>-10);
% if length(tmp)>1
%     s=norm(tmp,'fro')/length(VN);
% else
%     s=0;
% end
% if m-2*s>0%a>2*b
%     m
%     s
%     p=min(min(P(2:end,2:end)));
%     indAll=find(P<=p);
% else
%     p=P(end);
%     indAll=size(P,1)*size(P,2);
% end
% V=V(2:end,2:end);
k=Verify(V); % verify the rows
l=Verify(V'); % verify the columns
[m,n]=size(V);
% plot(V(2:m,n)-V(1:m-1,n));
% figure
% plot(V(m,2:n)-V(m-1,n-1));
V(end)
figure
plot(V(1:m,n))
figure
plot(V(m,1:n))
% std(V(1:m,n))
% std(V(m,1:n))
p=min(min(P(k,l)));
indAll=find(P<=p);


function k=Verify(VN)
VN(1,1:end)=0;
VN(1:end,1)=0;
[n,~]=size(VN);
a=mean(VN,2);
max(a(2:end));
% a(end)
c=a(2:n)-a(1:n-1);
% max(c)
% figure
% plot(c)
%ind=sort(c,'ascend');
% min(min(c))
m=mean(c);
s=std(c);
ind=find(c<-max(c));
if isempty(ind);
    k=n;
else
    k=find(c==min(c(ind)),1,'first')+1;
%     min(min(c))
end
% 
% rowMin=zeros(m,1);
% rowMax=zeros(m,1);
% for j=1:m
%     rowMin(j)=prctile(VN(j,:),5);
%     rowMax(j)=prctile(VN(j,:),95);
% end
% M=max(rowMax);

% 
% function [p,indAll]=MGCScaleVerify(V,P)
% % V=enlarge(V);
% % p=min(min(V(2:end,2:end)));
% % indAll=find(V<=p);
% % 
% % function V=enlarge(V)
% % [m,n]=size(V);
% % V(1:m,1)=1;
% % V(1,1:n)=1;
% % for i=2:m
% %     for j=2:n
% %         V(i,j)=V(i,j)*max((m-1)/(i-1),(n-1)/(j-1))^2;
% %     end
% % end
% 
% 
% % % An auxiliary function to verify and estimate the MGC optimal scale based
% % % on the p-values of all local correlations.
% VN=V(2:end,2:end);
% p1=Verify(VN); % verify the rows
% p2=Verify(VN'); % verify the columns
% p=min(p1,p2);
% p=min(p,V(end));
% indAll=find(V<=p);
% if isempty(indAll)
%     indAll=size(V,1)*size(V,2);
% end
% % p=V(indAll(end));
% 
% function p1=Verify(VN)
% [m,~]=size(VN);
% 
% rowTmp=median(VN,2)';
% for j=1:m
%     rowTmp(j)=prctile(VN(j,:),50);
% end
% %rowTmp=min(VN,[],2)';
% indK=find(rowTmp==min(rowTmp),1,'last');
% 
% % from the row with minimal median p-value, include adjacency rows
% % whose median p-value are significant
% tmpP=min(rowTmp)*m;
% % ss=indK;
% tt=indK;
% for i=indK-1:-1:1
%     tt=[i tt];
%     tmp=mean(rowTmp(tt));
%     if tmp*m/length(tt)<tmpP
%         tmpP=tmp*m/length(tt);
% %         ss=i;
%     else
%         break;
%     end
% end
% 
% % tmpP=min(rowTmp)*m;
% % se=indK;
% % tt=indK;
% for i=indK+1:m
%     tt=[tt i];
%     tmp=mean(rowTmp(tt));
%     if tmp*m/length(tt)<tmpP
%         tmpP=tmp*m/length(tt);
%         se=i;
%     else
%         break;
%     end
% end
% % 
% % indK=ss:1:se;
% indK=tt;
% p1=mean(rowTmp(indK))*m/length(indK);
% % check if the included rows satisfy a Bonferroni-type bound
% if p1>VN(end)
%     p1=rowTmp(end);
% end

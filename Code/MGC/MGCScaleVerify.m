function [p,indAll]=MGCScaleVerify(P,rep)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations
[p1,L]=Validate(P,rep);
[p2,K]=Validate(P',rep);
if p1<p2
    p=p1;
    P(:,L)=P(:,L)-1;
    indAll=find(P<=p-1);
else
    p=p2;
    P(K,:)=P(K,:)-1;
    indAll=find(P<=p-1);
end
% if (p<=0.05 && P(end)>0.05)
%     p
%     P(end)
% end

function [p,L]=Validate(P,rep)
n=size(P,2);
pL=median(P(2:end,1:end));
pL(1)=1;
thres2=max(0.0025,2/rep);

rs=ones(1,n);
for i=2:n
    tmp=diff(P(2:end,i));
   tmp(abs(tmp)<thres2)=0;
    tmp=tmp(1:n-3).*tmp(2:n-2);
    rs(i)=mean(tmp<0);
end
thres1=min(0.5-4*std(rs(2:end)),0.3);
pCount=zeros(n,1);
nn=max(1,round(0.025*n));
for i=2:n
    if rs(i)>thres1
        continue
    else
    tmp=(pL>pL(i)) | (rs>rs(i));
    tmp=find(tmp==1);
    tmp=[1;tmp';n+1];
    s1=find(tmp<i,1,'last');
    s2=find(tmp>i,1,'first');
    s1=tmp(s1)+1;
    s2=tmp(s2)-1;
       if s2>=min(i+nn,n) && s1<=max(2,i-nn)
           pCount(i)=(s2-s1+1)/(n-1);
       end
           pCount(i)=(pCount(i)>=max(0.1,5/(n-1)));
    end
end
L=find(pCount==1);
p=min(pL(L));
if isempty(L) || p>P(end)
    p=P(end);
    L=n;
else
    L=find(pL==p);
end


% thres1=10;
% thres2=0.5;
% [m,~]=size(P);
% coeR=ones(m,1);
% for i=2:m
%     d1=diff(P(i,2:end));
%     %d1(abs(d1)<max(2/rep,0.0025))=0;
%     d2=d1(1:length(d1)-1).*d1(2:length(d1));
%     coeR(i)=mean(d2<0);
% end
% % coeR=mean(P,2);
% %coeR=round(coeR*100)/100;
% %thres1=mean(coeR)-3*std(coeR)
% thres=max(coeR(end)-2*std(coeR(2:end)),thres1/100);
% %thres=0.1;
% 
% nn=max(ceil(length(coeR)*0.05),1);
% %nn=1;
% ind=find(coeR>thres);
% ind=[1;ind;m+1];
% t=diff(ind);
% t=find(t>2*nn+1);
% K=[];
% for i=1:length(t);
%     K=[K ind(t(i))+1:ind(t(i)+1)-1];
% end
% 
% if isempty(K)
%     K=m;
%     p=min(min(P(K,2:end),[],2));
% else
%     p=min(min(P(K,2:end),[],2));
% end
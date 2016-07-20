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
if (p<=0.05 && P(end)>0.05)
    p
    P(end)
%     figure
%     imagesc(P);
%     caxis([0 0.05])
end

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
% min(rs)

if min(pL)<min(2/rep,0.0025)
    tmp=find(pL>min(pL));
    tmp=max(diff(tmp));
    if tmp>2
        L=find(pL==min(pL));
        p=min(pL(L));
        return;
    end
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
           pCount(i)=(round(pCount(i)*100)/100>=max(0.15,5/(n-1)));
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
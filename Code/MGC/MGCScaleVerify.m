function [p,indAll]=MGCScaleVerify(P)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations
[p1,L]=Validate(P);
[p2,K]=Validate(P');
p=min(p1,p2);
% p=min(min(P(K,L)));
indAll=find(P<=p);
if p<=0.05 && P(end)>0.05
    p
    P(end)
end

function [p,L]=Validate(P)
n=size(P,1);
pL=median(P(2:end,2:end));
% pL=zeros(1,n-1);
% for i=2:n
%     pL(i-1)=prctile(P(:,i),25);
% end

if n<4
    p=max(pL);
    L=1:n;
    return;
end
nn=ceil(n*0.1);
%pL=[pL(3);pL';pL(n-3)];
pL=[zeros(nn,1);pL';zeros(nn,1)];
L=[];
p=1;
for i=2:n
    thres=pL(i);
    if sum(pL(i-1:i+2*nn-1)<=thres)==2*nn+1
        if thres<p;
           L=i;
           p=thres;
        end
    end
end
if isempty(L) || p>P(end)
    L=n;
    p=P(end);
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
function [p,indAll]=MGCScaleVerify(P)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations.
p1=verifyR(P);
p2=verifyR(P');
p=min(p1,p2);
indAll=find(P==p);
%  if p<=0.05 && P(end)>0.05
%     p
%     P(end)
% end

%p=max(max(P(indAll)));

function p=verifyR(P)
thres1=10;
thres2=0.3;
[m,~]=size(P);
coeR=ones(m,1);
for i=2:m
    d1=diff(P(i,2:end));
    d1(abs(d1)<=thres2/100)=0;
    d2=d1(1:length(d1)-1).*d1(2:length(d1));
    coeR(i)=mean(d2<0);
end
coeR=round(coeR*100)/100;

K=find(coeR<=thres1/100);
nn=max(floor(length(coeR)*thres2/10),1);
%nn=1;
k=[];
rk=round(prctile(coeR,(nn*2+1)/m*100)*100)/100;
coeR=[zeros(nn,1);coeR;zeros(nn,1)];
for i=1:length(K)
   if sum(coeR(K(i):K(i)+2*nn)<=rk)==2*nn+1 || sum(coeR(K(i):K(i)+2*nn)<=thres1/100)==2*nn+1
       k=[k K(i)];
   end
end
if isempty(k)
    p=P(end);
else
    p=max(min(P(k,2:end),[],2));
end
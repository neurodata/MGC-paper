function [p,indAll]=MGCScaleVerify(P)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations.
[coeR,coeC]=verifyAll(P);
coeR=[0;coeR;0]
coeC=[0;coeC;0]
thres=0.4;
% t1=(coeR<0.2);
% t2=(coeC<0.2);

ind=find(P<=min(min(P(2:end,2:end)))+0.01);
[m,n]=size(P);
[K,L]=ind2sub([m,n],ind);
indAll=[];
for i=1:length(K);
    k=K(i);
    l=L(i);
    %if sum(coeR(k:k+2)<thres)==3 || sum(coeC(l:l+2)<thres)==3
%     
     t1=coeR(k)+coeC(n);
     t2=coeR(m)+coeC(l);
     if t1<thres || t2<thres
        indAll=[indAll, sub2ind(size(P),k,l)];
    end
end
if isempty(indAll)
    indAll=(m)*(n);
end

p=max(max(P(indAll)));

function [coeR,coeC]=verifyAll(V)
[m,n]=size(V);
coeR=ones(m,1);
coeC=ones(n,1);
for i=2:m
    coeR(i)=verify(V(i,2:end));
end
for i=2:n
    coeC(i)=verify(V(2:end,i)');
end
% [~,~,coeR]=unique(coeR);
% [~,~,coeC]=unique(coeC);

function coe=verify(VC)
d1=diff(VC);
%d1(abs(d1)<0.005)=0;
d2=d1(1:length(d1)-1).*d1(2:length(d1));
coe=mean(d2<0);
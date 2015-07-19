function corrXY = rankDCorr(X,Y) %calculate rank dCorr
% Author: Cencheng Shen
% Implements the rank distance correlation from Shen, Jovo, CEP 2015?
n=size(X,1);
disA=disToRanks([X Y]); %change from distance to ranks of distance
X=disA(1:n,1:n);
Y=disA(1:n,n+1:2*n);
corrXY=zeros(n-1,1);
varX=zeros(n-1,1);
varY=zeros(n-1,1);

H=eye(n)-(1/n)*ones(n,n);
XH=X*H;
YH=Y*H;

for i=1:n
    for j=1:n
        tmp1=X(j,i);
        tmp2=Y(j,i);
        if tmp1~=0 && tmp2~=0
            corrXY(max(tmp1,tmp2):end)=corrXY(max(tmp1,tmp2):end)+XH(j,i)*YH(j,i);
            varX(tmp1:end)=varX(tmp1:end)+XH(j,i)*XH(j,i);
            varY(tmp2:end)=varY(tmp2:end)+YH(j,i)*YH(j,i);
        end
    end
end
corrXY=corrXY./sqrt(varX.*varY); %need to justify why we do not need negative value?

function disA=disToRanks(dis) %transform from distance to ranking, order from 0,...,n-1. Note that 1,...,n also works
n=size(dis,1);
[~,ind]=sort(dis(1:n,1:n));
[~,ind2]=sort(dis(1:n,n+1:2*n));
disA=zeros(n,2*n);
for i=1:n
    v=ind(2:n,i);
    w=ind2(2:n,i);
    for j=1:n-1
        disA(v(j),i)=j;
        disA(w(j),n+i)=j;
    end
end
function [mCorr,t] = mDistCorr(X,Y) %calculate modified dCorr
% Author: Cencheng Shen
% Implements the modified distance correlation from Szeley 2013
n=size(X,1);
H=eye(n)-(1/n)*ones(n,n);
mX=H*X*H;
mY=H*Y*H;
mX=n/(n-1)*(mX-X/n);
mY=n/(n-1)*(mY-Y/n);
for i=1:n
    mX(i,i)=0;
    mY(i,i)=0;
end
VXY=abs(sum(sum(mX.*(mY))));
VXX=abs(sum(sum(mX.*(mX))));
VYY=abs(sum(sum(mY.*(mY))));
meanX=sum(sum(X))/n^2;
meanY=sum(sum(Y))/n^2;
for i=1:n
    mX(i,i)=n/(n-1)*(mean(X(:,i))-meanX);
    mY(i,i)=n/(n-1)*(mean(Y(:,i))-meanY);
    VXY=VXY-2/(n-2)*mX(i,i)*mY(i,i);
    VXX=VXX-2/(n-2)*mX(i,i)*mX(i,i);
    VYY=VYY-2/(n-2)*mY(i,i)*mY(i,i);
end
VXY=VXY/n/(n-3);
VXX=VXX/n/(n-3);
VYY=VYY/n/(n-3);
if VXX<0 || VYY<0
    mCorr=0;
    t=0;
else
    mCorr=VXY/sqrt(VXX*VYY); %mCorr is asymptotically non-negative
    t=sqrt(n*(n-3)/2-1)*mCorr/sqrt(1-mCorr^2);
end
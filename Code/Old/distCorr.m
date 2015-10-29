function corr = distCorr(X,Y) %calculate dCorr
% Author: Cencheng Shen
% Implements the distance correlation from Szeley 2007
n=size(X,1);
H=eye(n)-(1/n)*ones(n,n);
X=H*X*H;
Y=H*Y*H;
corr=sum(sum(X.*(Y)));  %dist covariance is always non-negative
if corr<0
    corr=0;
else
    corr=sqrt(corr/(norm(X,'fro')*norm(Y,'fro')));
end
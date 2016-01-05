function corr = Mantel(X,Y) % Calculate Mantel statistic by Pearson's correlation
% Author: Cencheng Shen
n=size(X,1);
EX=sum(sum(X))/n/(n-1);
EY=sum(sum(Y))/n/(n-1);
X=X-EX;
Y=Y-EY;

corr=sum(sum(X.*Y))-sum(diag(X).*diag(Y));%-EX*EY;
var1=sum(sum(X.*X))-sum(diag(X).*diag(X));%-EX*EX;
var2=sum(sum(Y.*Y))-sum(diag(Y).*diag(Y));%-EY*EY;
corr=corr/sqrt(var1*var2);

function corr = Mantel(X,Y) % Calculate Mantel statistic by Pearson's correlation
% Author: Cencheng Shen
n=size(X,1);
EX=sum(sum(X))/n/(n-1);
EY=sum(sum(Y))/n/(n-1);
corr=sum(sum(X.*Y))/n/(n-1)-EX*EY;
var1=sum(sum(X.*X))/n/(n-1)-EX*EX;
var2=sum(sum(Y.*Y))/n/(n-1)-EY*EY;
corr=corr/sqrt(var1*var2); 
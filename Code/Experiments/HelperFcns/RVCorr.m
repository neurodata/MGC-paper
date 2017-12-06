function [corr] = RVCorr(x,y,option)
% Author: Cencheng Shen
% The main function that calculates all local correlation coefficients.
%
% The inputs are: 
% two symmetric distance matrices X and Y; or two n*d data matrices
% an option that specifies which global correlation to use, including 'mcor','dcor','mantel'.
%
% The outputs are all local correlations and all local variances.

if nargin<3
   option=0;
end
[n,d]=size(x);
x=x-repmat(mean(x,1),n,1);
y=y-repmat(mean(y,1),n,1);
cov=x'*y;
varx=x'*x;
vary=y'*y;
option=min(abs(option),d);
if option==0
   cov=trace(cov*cov');
   corr=cov/sqrt(trace(varx*varx)*trace(vary*vary));
else
   cov=sum(svds(cov,option).^2);
   corr=cov/sqrt(sum(svds(varx,option).^2)*sum(svds(vary,option).^2));
end

function [corr,varX,varY] = DCorr(A,B,option)
% Author: Cencheng Shen
% The main function that calculates all local correlation coefficients.
%
% The inputs are: 
% two symmetric distance matrices X and Y; or two n*d data matrices
% an option that specifies which global correlation to use, including 'mcor','dcor','mantel'.
%
% The outputs are all local correlations and all local variances.
if nargin < 3
    option='dcor'; % use dcor by default
end
if issymmetric(A)==false
    A=squareform(pdist(A));
end
if issymmetric(B)==false
    B=squareform(pdist(B));
end

[A,B]=MGCDistTransform(A,B,option,0);
[corr]=GlobalCov(A,B'); % compute all local corr / var statistics
[varX]=GlobalCov(A,A'); % compute all local corr / var statistics
[varY]=GlobalCov(B,B'); % compute all local corr / var statistics
% varX=diag(varX);
% varY=diag(varY);
corr=corr./real(sqrt(varX*varY'));

if varX<=0 || varY<=0
    corr=0;
end

if strcmp(option,'mantel')==true
    corr=abs(corr); % use absolute value for mantel coefficient
end

function [corr]=GlobalCov(A,B)
corr=sum(sum(A.*B));
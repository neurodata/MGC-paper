function [corr,varX,varY] = MGCLocalCorr(A,B,option)
% Author: Cencheng Shen
% The main function that calculates all local correlation coefficients.
%
% The inputs are: 
% two symmetric distance matrices X and Y; or two n*d data matrices
% an option that specifies which global correlation to use, including 'mcor','dcor','mantel'.
%
% The outputs are all local correlations and all local variances.
if nargin < 3
    option='mcor'; % use mcorr by default
end

if issymmetric(A)==false
    A=squareform(pdist(A));
end
if issymmetric(B)==false
    B=squareform(pdist(B));
end

[A,B,RX,RY]=MGCDistTransform(A,B,option);
[corr]=LocalCov(A,B,RX,RY); % compute all local corr / var statistics
[varX]=LocalCov(A,A',RX,RX'); % compute all local corr / var statistics
[varY]=LocalCov(B',B,RY',RY); % compute all local corr / var statistics
varX=diag(varX);
varY=diag(varY);
corr=corr./real(sqrt(varX*varY'));

for k=1:length(varX)
    if varX(k)<=0
        corr(k,:)=0;
    end
end
for l=1:length(varY)
    if varY(l)<=0
        corr(:,l)=0;
    end
end

function [corrXY]=LocalCov(A,B,RX,RY)
% An auxiliary function that computes all local correlations simultaneously in O(n^2)
[n,m]=size(A);nX=max(max(RX));nY=max(max(RY));
corrXY=zeros(nX,nY); %varX=zeros(1,nX); varY=zeros(1,nY);
EX=zeros(1,nX);EY=zeros(1,nY);

% summing up the entrywise product of A and B based on the ranks, which
% yields the local family of covariance and variances
for j=1:m
    for i=1:n
        a=A(i,j);b=B(i,j);k=RX(i,j);l=RY(i,j);
        corrXY(k,l)=corrXY(k,l)+a*b;
        EX(k)=EX(k)+a;
        EY(l)=EY(l)+b;
    end
end
for k=1:nX-1
    corrXY(k+1,1)=corrXY(k,1)+corrXY(k+1,1);
    EX(k+1)=EX(k)+EX(k+1);
end
for l=1:nY-1
    corrXY(1,l+1)=corrXY(1,l)+corrXY(1,l+1);
    EY(l+1)=EY(l)+EY(l+1);
end
for l=1:nY-1
    for k=1:nX-1
        corrXY(k+1,l+1)=corrXY(k+1,l)+corrXY(k,l+1)+corrXY(k+1,l+1)-corrXY(k,l);
    end
end

% normalize the covariance by the variances yields the local correlation
corrXY=(corrXY-EX'*EY/n^2);
corrXY(1,1:nY)=0;
corrXY(1:nX,1)=0;
% varX=varX-EX.^2/n^2; %%%new
% varY=varY-EY.^2/n^2; %%%new
% corrXY=corrXY./real(sqrt(varX'*varY));

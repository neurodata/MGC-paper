function [corr,varX,varY] = MGCLocalCorr(X,Y,option)
% Author: Cencheng Shen
% The main function that calculates all local correlation coefficients.
%
% The inputs are: 
% two distance matrices X and Y;
% an option that specifies which global correlation to use, including 'mcor','dcor','mantel'.
%
% The outputs are all local correlations and all local variances.
if nargin < 3
    option='mcor'; % use mcorr by default
end

[A,B,RX,RY]=MGCDistTransform(X,Y,option);
[corr,varX,varY]=LocalCorrelations(A,B,RX,RY); % compute all local corr / var statistics

function [corrXY,varX,varY]=LocalCorrelations(A,B,RX,RY)
% An auxiliary function that computes all local correlations simultaneously in O(n^2)
[n,m]=size(A);nX=max(max(RX));nY=max(max(RY));
corrXY=zeros(nX,nY); varX=zeros(1,nX); varY=zeros(1,nY);
EX=zeros(1,nX);EY=zeros(1,nY);

% summing up the entrywise product of A and B based on the ranks, which
% yields the local family of covariance and variances
for j=1:m
    for i=1:n
        a=A(i,j);b=B(i,j);k=RX(i,j);l=RY(i,j);
        corrXY(k,l)=corrXY(k,l)+a*b;
        varX(k)=varX(k)+a^2;
        varY(l)=varY(l)+b^2;
        EX(k)=EX(k)+a;
        EY(l)=EY(l)+b;
    end
end
for k=1:nX-1
    corrXY(k+1,1)=corrXY(k,1)+corrXY(k+1,1);
    varX(k+1)=varX(k)+varX(k+1);
    EX(k+1)=EX(k)+EX(k+1);
end
for l=1:nY-1
    corrXY(1,l+1)=corrXY(1,l)+corrXY(1,l+1);
    varY(l+1)=varY(l)+varY(l+1);
    EY(l+1)=EY(l)+EY(l+1);
end
for l=1:nY-1
    for k=1:nX-1
        corrXY(k+1,l+1)=corrXY(k+1,l)+corrXY(k,l+1)+corrXY(k+1,l+1)-corrXY(k,l);
    end
end

% normalize the covariance by the variances yields the local correlation
corrXY=(corrXY-EX'*EY/n^2);
varX=varX-EX.^2/n^2;
varY=varY-EY.^2/n^2;
corrXY=corrXY./real(sqrt(varX'*varY));

% set any local correlation to 0 if any corresponding local variance is no larger than 0
for k=1:nX
    if varX(k)<=0
        corrXY(k,:)=0;
    end
end
for l=1:nY
    if varY(l)<=0
        corrXY(:,l)=0;
    end
end

% function [corrXY,varX,varY]=GlobalCorrelation(A,B)
% corrXY=sum(sum(A.*B));
% varX=sum(sum(A.*A));
% varY=sum(sum(B.*B));
%
% % Normalizing the covariance by the variances yields the local correlation.
% if varX<=0 || varY<=0
%     corrXY=0;
% else
%     corrXY=corrXY/real(sqrt(varX*varY'));
% end

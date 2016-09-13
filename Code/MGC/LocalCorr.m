function [corrXY,varX,varY, weightXY] = LocalCorr(X,Y,option, optimalInd)
% Author: Cencheng Shen
% The main function that calculates all local correlation coefficients.
%
% The inputs are: 
% two distance matrices X and Y;
% an option that specifies which global correlation to use, including 'mcor','dcor','mantel'.
%
% The outputs are all local correlations and all local variances.
%
% Alternatively, specifying optimalInd by a matrix single index will return a
% weightXY matrix that shows the contribution of each distance entries to
% the eventual local distance correlation at the given index.
if nargin < 3
    option='mcor'; % use mcorr by default
end
if nargin < 4
    optimalInd=[];
end
n=size(X,1);
disRank=[DistRanks(X) DistRanks(Y)]; % sort distances within columns

% depending on the choice of the global correlation, properly center the distance matrices
A=DistCentering(X,option);
B=DistCentering(Y,option);

RX=disRank(1:n,1:n); % the column ranks for X
RY=disRank(1:n,n+1:2*n); % the column ranks for Y
if isempty(optimalInd)
    [corrXY,varX,varY]=LocalCorrelationComputation(A,B',RX,RY'); % compute all local corr / var statistics
    weightXY=[];
else
    corrXY=[]; varX=[]; varY=[];
    weightXY=LocalWeightComputation(A,B',RX,RY', optimalInd); % compute distance entry contributions
end

function [corrXY,varX,varY]=LocalCorrelationComputation(A,B,RX,RY)
% An auxiliary function that computes all local correlations simultaneously in O(n^2)
n=size(A,1);nX=max(max(RX));nY=max(max(RY));
corrXY=zeros(nX,nY); varX=zeros(1,nX); varY=zeros(1,nY);
EX=zeros(1,nX);EY=zeros(1,nY);

% summing up the entrywise product of A and B based on the ranks, which
% yields the local family of covariance and variances
for j=1:n
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
corrXY=(corrXY)./real(sqrt(varX'*varY));

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

function weightXY=LocalWeightComputation(A,B,RX,RY,ind)
% An auxiliary function that computes the contributions of each distance entries to
% the local distance correlation at a given scale.
nX=max(max(RX));nY=max(max(RY));
[k,l]=ind2sub([nX,nY],ind);
RX=(RX>k);
RY=(RY>l);
A(RX)=0;
B(RY)=0;
weightXY=(A-mean(mean(A))).*(B-mean(mean(B)));

% function [corrXY,varX,varY]=GlobalComputation(A,B)
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

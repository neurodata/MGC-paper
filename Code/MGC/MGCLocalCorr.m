%% Compute all local correlation coefficients in O(n^2 log n)
%%
%% @param X is a distance matrix or a n*d data matrix; if it is not a square matrix with zeros on diagonal, it is treated as n*d data; 
%% @param Y is a second distance matrix or a n*d data matrix, with the same distance matrix check as X;
%% @param option is a string that specifies which global correlation to build up-on, including 'mgc','dcor','mantel', and 'rank'.
%%
%% @return A list contains the following output:
%% @return corr consists of all local correlations within [-1,1] by double matrix index;
%% @return varX contains all local variances for X; varY contains all local covariances for Y.
%%
%% @export
%% 
function [corr,varX,varY] = MGCLocalCorr(X,Y,option)
if nargin < 3
    option='mgc'; % use mgc by default
end

% Use the data size and diagonal element to determine if the given data is a distance matrix or not
if size(X,1)~=size(X,2) || sum(diag(X).^2)>0
    X=squareform(pdist(X));
%     disp('The first data is not a Euclidean distance matrix; transformed to distance matrix instead.')
end
if size(Y,1)~=size(Y,2) || sum(diag(Y).^2)>0
    Y=squareform(pdist(Y));
%     disp('The second data is not a Euclidean distance matrix; transformed to distance matrix instead.')
end

[A,B,RX,RY]=MGCDistTransform(X,Y,option);
[corr]=LocalCov(A,B',RX,RY'); % compute all local covariances
[varX]=LocalCov(A,A',RX,RX'); % compute local variances for first data
[varY]=LocalCov(B,B',RY,RY'); % compute local variances for second data
varX=diag(varX);
varY=diag(varY);
corr=corr./real(sqrt(varX*varY'));
corr(corr>1)=1; % avoid computational issue that may cause a few local corr to be negligably larger than 1

% set any local correlation to 0 if any corresponding local variance is no larger than 0
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

%% An auxiliary function that computes all local correlations simultaneously in O(n^2).
%%
%% @param A is a properly transformed distance matrix;
%% @param B is the second distance matrix properly transformed;
%% @param RX is the column-ranking matrix of A;
%% @param RY is the column-ranking matrix of B.
%%
%% @return covXY is all local covariances computed iteratively.
%% 
function [covXY]=LocalCov(A,B,RX,RY)
n=size(A,1);nX=max(max(RX));nY=max(max(RY));
covXY=zeros(nX,nY); %varX=zeros(1,nX); varY=zeros(1,nY);
EX=zeros(1,nX);EY=zeros(1,nY);

% summing up the entrywise product of A and B based on the ranks, which
% yields the local family of covariance and variances
for j=1:n
    for i=1:n
        a=A(i,j);b=B(i,j);k=RX(i,j);l=RY(i,j);
        covXY(k,l)=covXY(k,l)+a*b;
        EX(k)=EX(k)+a;
        EY(l)=EY(l)+b;
    end
end
for k=1:nX-1
    covXY(k+1,1)=covXY(k,1)+covXY(k+1,1);
    EX(k+1)=EX(k)+EX(k+1);
end
for l=1:nY-1
    covXY(1,l+1)=covXY(1,l)+covXY(1,l+1);
    EY(l+1)=EY(l)+EY(l+1);
end
for l=1:nY-1
    for k=1:nX-1
        covXY(k+1,l+1)=covXY(k+1,l)+covXY(k,l+1)+covXY(k+1,l+1)-covXY(k,l);
    end
end

% normalize the covariance by the variances yields the local correlation
covXY=(covXY-EX'*EY/n/(n));
% covXY(1,1:nY)=0; % local cov without any neighbor is meaningless and set to 0 instead
% covXY(1:nX,1)=0;
% varX=varX-EX.^2/n^2; %%%new
% varY=varY-EY.^2/n^2; %%%new
% corrXY=corrXY./real(sqrt(varX'*varY));

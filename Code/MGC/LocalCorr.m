function [corrXY,varX,varY] = LocalCorr(X,Y,option)
% Author: Cencheng Shen
% The main function that calculates all local correlation coefficients.
%
% The inputs are: 
% two distance matrices X and Y;
% an option that specifies which global correlation to use, set to 1,2,3 for dcorr / mcorr / Mantel.
%
% The outputs are all local correlations and all local variances.
if nargin < 3
    option=2; % use mcorr by default
end
n=size(X,1);
disRank=[disToRanks(X) disToRanks(Y)]; % sort distances within columns

% depending on the choice of the global correlation, properly center the
% distance matrices
A=Centering(X,option);
B=Centering(Y,option);

RX=disRank(1:n,1:n); % the ranks for X
RY=disRank(1:n,n+1:2*n); % the ranks for Y
[corrXY,varX,varY]=LocalComputation(A',B,RX,RY); % compute all local corr / var statistics

function [A]=Centering(X,option)
% An auxiliary function that properly centers the distance matrix X,
% depending on the choice of global corr.
n=size(X,1);
if option==3
    % centering for Mantel
    EX=sum(sum(X))/n/(n-1);
    A=X-EX;
    for j=1:n
        A(j,j)=0;
    end
else
    % centering for dcorr / mcorr, which uses single centering rather than
    % the original double centering
    A=X-repmat(mean(X,1),n,1);
    %A=X-repmat(mean(X,1),n,1)-repmat(mean(X,2),1,n)+mean(mean(X)); % this is original double-centering
    
    % for mcorr, further adjust the centered matrices to remove high-dimensional bias
    if option==2
        A=A-X/n;
        %%%   meanX=sum(sum(X))/n^2;
        for j=1:n
            A(j,j)=0;
            %%% A(j,j)=sqrt(2/(n-2))*(mean(X(:,j))-meanX)*1i;  % the original diagonal modification of mcorr
        end
    end
end

function [corrXY,varX,varY]=LocalComputation(A,B,RX,RY)
% An auxiliary function that computes all local correlations simultaneously
% in O(n^2)
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

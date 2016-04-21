function [corrXY,varX,varY] = LocalGraphCorr(X,Y,option,disRank) % Calculate local variants of mcorr/dcorr/Mantel
% Author: Cencheng Shen
% Implements local graph correlation from Shen, Jovo, CEP 2016.
%
% By specifying option=1, 2, or 3, it calculates the local correlations of dcorr, mcorr, and Mantel
% Specifying the rank matrix by a size n * 2n disRank can save the sorting.
if nargin < 3
    option=1; % By default use dcorr
end
if nargin < 4
    disRank=[disToRanks(X) disToRanks(Y)]; % Sort distances within columns, if the ranks are not provided
end
n=size(disRank,1);
RX=disRank(1:n,1:n); % The ranks for X
RY=disRank(1:n,n+1:2*n); % The ranks for Y

% The output contains the local variants of correlations and variances
corrXY=zeros(n,n); varX=zeros(n,1); varY=zeros(n,1);

% Depending on the choice of the global test, calculate the entries of A and B
% accordingly for late multiplication.
if option~=3
    % Double centering for dcorr/mcorr
    A=doubleCentering(X);
    B=doubleCentering(Y);
    % For mcorr, further adjust the double centered matrices to remove high-dimensional bias
    if option==2
        A=A-X/n;
        B=B-Y/n;
%         meanX=sum(sum(X))/n^2;
%         meanY=sum(sum(Y))/n^2;
        % The diagonals of mcorr are set to zero, instead of the original formulation of mcorr
        for j=1:n
            A(j,j)=0;
            B(j,j)=0;
% %             The original diagonal modification of mcorr
%             A(j,j)=sqrt(2/(n-2))*(mean(X(:,j))-meanX)*1i;
%             B(j,j)=sqrt(2/(n-2))*(mean(Y(:,j))-meanY)*1i;
        end
    end
else
    % Single centering for Mantel
    EX=sum(sum(X))/n/(n-1);
    EY=sum(sum(Y))/n/(n-1);
    A=X-EX;
    B=Y-EY;
    % Mantel does not use diagonal entries, which is equivalent to set them zero
    for j=1:n
        A(j,j)=0;
        B(j,j)=0;
    end
end

% Summing up the entrywise product of A and B based on the ranks, which
% yields the local family of covariance and variances
for j=1:n
    for i=1:n
        a=A(i,j);
        b=B(i,j);
        % If there are ties, set all rank 0 entries to the diagonal entry
        if (RX(i,j)==0)
            a=A(j,j);
        end
        if (RY(i,j)==0)
            b=B(j,j);
        end
        tmp1=RX(i,j)+1;
        tmp2=RY(i,j)+1;
        corrXY(tmp1,tmp2)=corrXY(tmp1,tmp2)+a*b;
        varX(tmp1)=varX(tmp1)+a^2;
        varY(tmp2)=varY(tmp2)+b^2;
    end
end
for j=1:n-1
    corrXY(1,j+1)=corrXY(1,j)+corrXY(1,j+1);
    corrXY(j+1,1)=corrXY(j,1)+corrXY(j+1,1);
    varX(j+1)=varX(j)+varX(j+1);
    varY(j+1)=varY(j)+varY(j+1);
end
for j=1:n-1
    for i=1:n-1
        corrXY(i+1,j+1)=corrXY(i+1,j)+corrXY(i,j+1)+corrXY(i+1,j+1)-corrXY(i,j);
    end
end

% Normalizing the covariance by the variances yields the local correlation.
corrXY=corrXY./real(sqrt(varX*varY'));

% Set any local correlation to 0 if any corresponding local variance is no larger than 0
for j=1:length(varX)
    if varX(j)<=0
        corrXY(j,:)=0;
    end
    if varY(j)<=0
        corrXY(:,j)=0;
    end
end

function A=doubleCentering(X)
%H=eye(n)-(1/n)*ones(n,n);
%A=H*X*H;
n=size(X,1);
A=mean(X,1);
A=repmat(A,n,1)+repmat(mean(X,2),1,n)-mean(A);
%A=repmat(A,n,1)+repmat(A',1,n)-mean(A);
A=X-A;

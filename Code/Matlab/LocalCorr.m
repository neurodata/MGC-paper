function [corrXY,varX,varY] = LocalCorr(X,Y,option,disRank) % Calculate local variants of mcorr/dcorr/Mantel
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
n=size(X,1);

% Depending on the choice of the global test, calculate the entries of A and B
% accordingly for late multiplication.
A=Centering(X,option);
B=Centering(Y,option);

%   [corrXY,varX,varY]=GlobalComputation(A',B);
RX=disRank(1:n,1:n); % The ranks for X
RY=disRank(1:n,n+1:2*n); % The ranks for Y
[corrXY,varX,varY]=LocalComputation(A',B,RX,RY);

function [A]=Centering(X,option)
n=size(X,1);
if option==3
    EX=sum(sum(X))/n/(n-1);
    A=X-EX;
    % Mantel does not use diagonal entries, which is equivalent to set them zero
    for j=1:n
        A(j,j)=0;
    end
else
    % Double centering for dcorr/mcorr
    %A=X-repmat(mean(X,1),n,1)-repmat(mean(X,2),1,n)+mean(mean(A));
    A=X-repmat(mean(X,1),n,1);
    % For mcorr, further adjust the double centered matrices to remove high-dimensional bias
    if option==2
        A=A-X/n;
        %         meanX=sum(sum(X))/n^2;
        % The diagonals of mcorr are set to zero, instead of the original formulation of mcorr
        for j=1:n
            A(j,j)=0;
            % %             The original diagonal modification of mcorr
            %             A(j,j)=sqrt(2/(n-2))*(mean(X(:,j))-meanX)*1i;
        end
    end
end

function [corrXY,varX,varY]=LocalComputation(A,B,RX,RY)
% The output contains the local variants of correlations and variances
n=size(A,1);
nX=max(max(RX));nY=max(max(RY));
corrXY=zeros(nX,nY); varX=zeros(nX,1); varY=zeros(nY,1);
EX=zeros(nX,1);
EY=zeros(nY,1);

% Summing up the entrywise product of A and B based on the ranks, which
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
% Normalizing the covariance by the variances yields the local correlation.
corrXY=(corrXY-EX*EY'/n^2);%;
varX=varX-EX.*EX/n^2;
varY=varY-EY.*EY/n^2;
corrXY=(corrXY)./real(sqrt(varX*varY'));

% Set any local correlation to 0 if any corresponding local variance is no larger than 0
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

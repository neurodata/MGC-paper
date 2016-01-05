function [corrXY,varX,varY] = LocalGraphCorr(X,Y,option, disRank) % Calculate local dCorr
% Author: Cencheng Shen
% Implements local distance correlation from Shen, Jovo, CEP 2015.
% By specifying optionModified=1, 2, or 3, it calculates the LGC statistic
% with respect to mcorr, dcorr, and Mantel respectively
if nargin < 3
    option=1; % By default use mcorr
end
if nargin < 4
    disRank=[disToRanks(X) disToRanks(Y)]; % Sort distances within columns, if the ranks are not provided
end
n=size(disRank,1);
RX=disRank(1:n,1:n);
RY=disRank(1:n,n+1:2*n);
corrXY=zeros(n,n);
varX=zeros(n,1);
varY=zeros(n,1);

% Double centering of distance matrices
if option~=3
    H=eye(n)-(1/n)*ones(n,n);
    A=H*X*H;
    B=H*Y*H;
    % Adjust the centering for modified dCorr
    if option==1
        A=A-X/n;
        B=B-Y/n;
        %                 meanX=sum(sum(X))/n^2;
        %                 meanY=sum(sum(Y))/n^2;
        for j=1:n
            %                         A(j,j)=sqrt(2/(n-2))*(mean(X(:,j))-meanX)*1i;
            %                         B(j,j)=sqrt(2/(n-2))*(mean(Y(:,j))-meanY)*1i;
            A(j,j)=0;
            B(j,j)=0;
        end
    end
else
    EX=sum(sum(X))/n/(n-1);
    EY=sum(sum(Y))/n/(n-1);
    A=X-EX;
    B=Y-EY;
    for j=1:n
        A(j,j)=0;
        B(j,j)=0;
    end
end

% Summing up the product of distances, which yields dCov and dVar
for j=1:n
    for i=1:n
        a=A(i,j);
        b=B(i,j);
        if (RX(i,j)==0)
            a=A(j,j);
        end
        if (RY(i,j)==0)
            b=B(j,j);
        end
        tmp1=RX(i,j)+1;
        tmp2=RY(i,j)+1;
        corrXY(tmp1:end, tmp2:end)=corrXY(tmp1:end, tmp2:end)+a*b;
        varX(tmp1:end)=varX(tmp1:end)+a^2;
        varY(tmp2:end)=varY(tmp2:end)+b^2;
    end
end

% Normalizing dCov by dVar yields dCorr
corrXY=corrXY./real(sqrt(varX*varY'));

% Set dCorr to 0 if any dVar is no larger than 0
for j=1:length(varX)
    if varX(j)<=0
        corrXY(j,:)=0;
    end
    if varY(j)<=0
        corrXY(:,j)=0;
    end
end
% % The original dCorr is defined as the square root of the previous calculated dCorr; but square root or not does not affect testing at all.
% if option==2
%     corrXY=real(sqrt(corrXY));
% end
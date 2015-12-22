function [corrXY,varX,varY] = localDCorr(X,Y,optionModified, disRank) % Calculate local dCorr
% Author: Cencheng Shen
% Implements local distance correlation from Shen, Jovo, CEP 2015.
% By specifying optionModified=0 or 1, it incorporates knn into original dcorr and modified dcorr
if nargin < 3
    optionModified=0; % By default use original dCorr
end
if nargin<4
    disRank=[disToRanks(X) disToRanks(Y)]; % Sort distances within columns, if the ranks are not provided
end
n=size(disRank,1);
RX=disRank(1:n,1:n);
RY=disRank(1:n,n+1:2*n);
corrXY=zeros(n,n);
varX=zeros(n,1);
varY=zeros(n,1);

% Double centering of distance matrices
H=eye(n)-(1/n)*ones(n,n);
XH=H*X*H;
YH=H*Y*H;
% Adjust the centering for modified dCorr
if optionModified==(1)
    XH=n/(n-1)*(XH-X/n);
    YH=n/(n-1)*(YH-Y/n);
    for i=1:n
        XH(i,i)=0;
        YH(i,i)=0;
    end
end

% Summing up the product of distances, which yields dCov and dVar
for i=1:n
    if optionModified==0
        corrXY=corrXY+XH(i,i)*YH(i,i);
        varX=varX+XH(i,i)*XH(i,i);
        varY=varY+YH(i,i)*YH(i,i);
    end
    for j=1:n
        if RX(j,i)>0&&RY(j,i)>0%i~=j
            tmp1=RX(j,i)+1;
            tmp2=RY(j,i)+1;
            corrXY(tmp1:end, tmp2:end)=corrXY(tmp1:end, tmp2:end)+XH(j,i)*YH(j,i);
            varX(tmp1:end)=varX(tmp1:end)+XH(j,i)*XH(j,i);
            varY(tmp2:end)=varY(tmp2:end)+YH(j,i)*YH(j,i);
        end
    end
end

% Further adjust the diagonal products of modified dCov
if optionModified==1
    meanX=sum(sum(X))/n^2;
    meanY=sum(sum(Y))/n^2;
    for i=1:n
        XH(i,i)=n/(n-1)*(mean(X(:,i))-meanX);
        YH(i,i)=n/(n-1)*(mean(Y(:,i))-meanY);
        corrXY=corrXY-2/(n-2)*XH(i,i)*YH(i,i);
        varX=varX-2/(n-2)*XH(i,i)*XH(i,i);
        varY=varY-2/(n-2)*YH(i,i)*YH(i,i);
    end
    corrXY=corrXY/n/(n-3);
    varX=varX/n/(n-3);
    varY=varY/n/(n-3);
end
% Normalizing dCov by dVar yields dCorr
corrXY=corrXY./real(sqrt(varX*varY'));

% Set dCorr to 0 if any dVar is no larger than 0
for i=1:n-1
    if varX(i)<=0
        corrXY(i,:)=0;
    end
    if varY(i)<=0
        corrXY(:,i)=0;
    end
end
% The original dCorr is defined as the square root of the previous calculated dCorr; but square root or not does not affect testing at all.
if optionModified==0
    corrXY=real(sqrt(corrXY));
end
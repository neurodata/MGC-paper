function [corrXY,varX,varY] = rankDCorr(X,Y,disRank,option) %calculate rank dCorr
% Author: Cencheng Shen
% Implements the rank distance correlation from Shen, Jovo, CEP 2015.
% By specifying option=1,2,3, it incorporates knn into rank dcorr, dcorr,
% and modified dcorr
if nargin < 3
    disRank=[disToRanks(X) disToRanks(Y)]; %change from distance to ranks of distance
end
if nargin<4
    option=1;
end
n=size(disRank,1);

RX=disRank(1:n,1:n);
RY=disRank(1:n,n+1:2*n);
corrXY=zeros(n,1);
varX=zeros(n,1);
varY=zeros(n,1);

H=eye(n)-(1/n)*ones(n,n);
XH=H*X*H;
YH=H*Y*H;

if option==(1)
    H2=(1/(n))*ones(n,n);
    XH=XH-H2*X*H2;
    YH=YH-H2*Y*H2;
    XH=n/(n-2)*(XH-2*X/n);
    YH=n/(n-2)*(YH-2*Y/n);
    for i=1:n
        XH(i,i)=0;
        YH(i,i)=0;
    end
end
if option==(3)
    XH=n/(n-1)*(XH-X/n);
    YH=n/(n-1)*(YH-Y/n);
    for i=1:n
        XH(i,i)=0;
        YH(i,i)=0;
    end
end

for i=1:n
    if option==2 %|| option==1
        corrXY=corrXY+XH(i,i)*YH(i,i);
        varX=varX+XH(i,i)*XH(i,i);
        varY=varY+YH(i,i)*YH(i,i);
    end
    for j=1:n
        if i~=j
            tmp1=RX(j,i)+1;
            tmp2=RY(j,i)+1;
            tmpMax=max(tmp1,tmp2);
            corrXY(tmpMax:end)=corrXY(tmpMax:end)+XH(j,i)*YH(j,i);
            varX(tmp1:end)=varX(tmp1:end)+XH(j,i)*XH(j,i);
            varY(tmp2:end)=varY(tmp2:end)+YH(j,i)*YH(j,i);
        end
    end
end

if option==3 %|| option==1
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
    corrXY=corrXY./sqrt(varX.*varY); %mCorr is asymptotically non-negative
else
    corrXY=corrXY./sqrt(varX.*varY);
    if option==2
        corrXY=sqrt(corrXY);
    end
end

for i=1:n-1
    if varX(i)<=0 || varY(i)<=0
        corrXY(i)=0;
    end
end
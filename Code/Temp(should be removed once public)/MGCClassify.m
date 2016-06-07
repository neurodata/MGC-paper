function [LabelTsn]=MGCClassify(dist,TrainInd, TestInd, LabelTrn)

% n=size(dist,1);
tsn=length(TestInd);
trn=length(TrainInd);
n=size(dist,1);
rep=500;
option=0;
if length(LabelTrn)~=trn || trn+tsn>n
    disp('Incorrect Size for Inputs');
end
LabelTsn=zeros(tsn,1);

[~,rowMean,allMean]=doubleCentering(dist);
distTsn=dist(:,TestInd);
distTsn=distTsn-repmat(rowMean',1,tsn)-repmat(mean(distTsn,1),n,1)+allMean;
distTsn=distTsn(TrainInd,:);
distTrn=dist(TrainInd,TrainInd);

% distTrn=dist(TrainInd,TrainInd);
% [~,rowMean,allMean]=doubleCentering(distTrn);
% distTsn=dist(TrainInd,TestInd);
% distTsn=distTsn-repmat(rowMean',1,tsn)-repmat(mean(distTsn,1),trn,1)+allMean;

[LabelCount,pos]=unique(LabelTrn);
distLabel=squareform(pdist(LabelTrn));
distLabel(distLabel>0)=1;
distLabel=doubleCentering(distLabel);
testLabel=distLabel(:,pos);

if option~=0
    [~,p1All,~,testAll,ind]=MGCPermutationTest(distTrn,distLabel,rep,option);
    [k,~]=ind2sub(size(p1All),ind);
else
    k=20;
end
% [~,varX,varY] = LocalGraphCorr(distTrn,distLabel,option);% Calculate local variants of mcorr/dcorr/Mantel
% varX=varX(k);varY=varY(l);

for i=1:tsn
    tmpTsn=distTsn(:,i);
    rkTsn=disToRanks(tmpTsn);
    tmpTsn(rkTsn>=k)=0;
    
    tmpM=-1000;
    tmpLabel=LabelCount(1);
    for j=1:size(testLabel,2)
        mm=sum(tmpTsn.*testLabel(:,j));
        if (mm>tmpM)
            tmpM=mm;
            tmpLabel=LabelCount(j);
        end
    end
    LabelTsn(i)=tmpLabel;
end

function [A, rowMean, allMean]=doubleCentering(X)
%H=eye(n)-(1/n)*ones(n,n);
%A=H*X*H;
n=size(X,1);
rowMean=mean(X,1);
allMean=mean(rowMean);
A=repmat(rowMean,n,1)+repmat(mean(X,2),1,n)-allMean;
%A=repmat(A,n,1)+repmat(A',1,n)-mean(A);
A=X-A;


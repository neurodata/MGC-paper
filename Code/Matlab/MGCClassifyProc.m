function [error]=MGCClassifyProc(Data,Label,trn,tsn)

% [error1]=MGCClassifyProc(Adj,Label);
% n=size(Adj,1);
% SRC_Euc_Test(Adj,ceil(n/2),floor(n/2),50,'resultBeta2SBMC',Label)
% [error1]=MGCClassifyProc(fea(per(1,:),:),gnd(per(1,:)));
[n,~]=size(Data);
if nargin<3 || trn+tsn>n;
    trn=ceil(n/2);
    tsn=floor(n/2);
end
dist=squareform(pdist(Data));

TrainInd=n-trn+1:n;
% distTrn=dist(TrainInd,TrainInd);
LabelTrn=Label(TrainInd);
TestInd=1:tsn;
% distTsn=dist(TrainInd,TestInd);
LabelTsn=Label(TestInd);

LabelTst=MGCClassify(dist,TrainInd, TestInd, LabelTrn);
error=mean(LabelTst~=LabelTsn);
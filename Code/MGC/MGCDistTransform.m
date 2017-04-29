function [A,B,RX,RY]=MGCDistTransform(X,Y,option)

if nargin<3
    option='mcor';
end
% depending on the choice of the global correlation, properly center the distance matrices
A=DistCentering(X,option);
B=DistCentering(Y,option);

RX=DistRanks(X); % the column ranks for X
RY=DistRanks(Y); % the column ranks for Y

if strcmp(option,'mcor') || strcmp(option,'dcor')
    B=B';
    RY=RY';
end

function [disRank]=DistRanks(dis)
% An auxiliary function that sorts the distance entries within each column by ascending order.
%
% The input is assumed to be a distance matrix.
%
% The output is column-wise rank, ordered from 1,...,n.
%
% For ties, the minimum ranking is used,
% e.g. if there are repeating distance entries, the order is like 1,2,3,3,4,..,n-1.

[n,m]=size(dis);
disRank=zeros(n,m);
for i=1:m
    [~,~,a]=unique(dis(:,i));
    disRank(:,i)=a;
end

function [A]=DistCentering(X,option)
% An auxiliary function that properly centers the distance matrix X,
% depending on the choice of global corr.
[n,m]=size(X);

% centering for distance correlation / modified distance correlation /
% Mantel coefficient
switch option
    case 'dcor' % single centering of dcor
        EX=repmat(mean(X,1),n,1); % column centering
    case 'mcor' % single centering of mcor
        EX=repmat(sum(X,1)/n,n,1);
        EX=EX+X/n;
    case 'mantel'
        EX=sum(sum(X))/n/(n-1);
    case 'dcorDouble' % original double centering of dcor
        EX=repmat(mean(X,1),n,1)+repmat(mean(X,2),1,n)-mean(mean(X));
    case 'mcorDouble' % original double centering of mcor
        EX=repmat(mean(X,1),n,1)+repmat(mean(X,2),1,n)-mean(mean(X));
        EX=EX+X/n;
end
A=X-EX;

% for mcor or Mantel, exclude the diagonal entries.
% This is a simpler diagonal modification than the original mcor
if strcmp(option,'dcor')==false && strcmp(option,'dcorDouble')==false
    %%% meanX=sum(sum(X))/n^2;
    for j=1:m
        A(j,j)=0;
        %%% A(j,j)=sqrt(2/(n-2))*(mean(X(:,j))-meanX)*1i;  % the original diagonal modification of mcorr
    end
end
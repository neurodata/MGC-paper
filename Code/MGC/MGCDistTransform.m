%% Transform the distance matrices, with column-wise ranking if needed.
%%
%% @param X is a distance matrix;
%% @param Y is a second distance matrix;
%% @param option is a string that specifies which global correlation to build up-on, including 'mgc','dcor','mantel', and 'rank';
%% @param optionRk is a string that specifies whether ranking within column is computed or not.
%%
%% @return A list contains the following output:
%% @return A and B are the centered distance matrices;
%% @return RX and RY are the column rank matrices of X and Y respectively.
%%
%% @export
%% 
function [A,B,RX,RY]=MGCDistTransform(X,Y,option,optionRk)
if nargin<3
    option='mgc'; % default to mgc transform
end
if nargin<4 || strcmp(option,'rank')
    optionRk=1; % do ranking or not, 0 to no ranking
end

% Depending on the choice of the global correlation, properly transform
% each distance matrix
if strcmp(option,'ReducedRank')
    RX=DistRanks(X);
    RY=DistRanks(Y);
    A=RX;
    B=RY;
    n=size(X,1);
    for i=1:n
        [~,aInd]=sort(A(:,i));
        b=sort(B(:,i));
        A(aInd,i)=b;
    end
    option='mgc'; % default to mgc transform
    [A,~]=DistCentering(A,option,0);
    [B,~]=DistCentering(B,option,0);
else
    [A,RX]=DistCentering(X,option,optionRk);
    [B,RY]=DistCentering(Y,option,optionRk);
end

%% An auxiliary function that properly transforms the distance matrix X
%%
%% @param X is a symmetric distance matrix;
%% @param option is a string that specifies which global correlation to build up-on, including 'mgc','dcor','mantel', and 'rank';
%% @param optionRk is a string that specifies whether ranking within column is computed or not.
%%
%% @return A list contains the following output:
%% @return A is the centered distance matrices;
%% @return RX is the column rank matrices of X.
%%
function [A,RX]=DistCentering(X,option,optionRk)
[n]=size(X,1);
if optionRk~=0
    RX=DistRanks(X); % the column ranks for X
else
    RX=zeros(n,n);
end

% If rank transformation, take X as the rank matrix then do default mgc
% transform
if strcmp(option,'rank')
    X=RX;
end

switch option
    case 'dcor' % unbiased dcor transform
        EX=repmat(sum(X,1)/(n-2),n,1)+repmat(sum(X,2)/(n-2),1,n)-sum(sum(X))/(n-1)/(n-2);
    case 'mantel' % mantel transform
        EX=sum(sum(X))/n/(n-1);
    otherwise % Default mgc transform
        EX=repmat(sum(X,1)/(n-1),n,1);
        %     case 'dcor' % original dcor that is biased
        %         EX=repmat(mean(X,1),n,1)+repmat(mean(X,2),1,n)-mean(mean(X));
        %     case 'dcor' % single centering of dcor
        %         EX=repmat(mean(X,1),n,1); % column centering
        %     case 'mcor' % single centering of mcor
        %         EX=repmat(sum(X,1)/n,n,1);
        %         EX=EX+X/n;
end
A=X-EX;

% The diagonal entries are always excluded
for j=1:n
    A(j,j)=0;
end

%% An auxiliary function that sorts the entries within each column by ascending order:
%% For ties, the minimum ranking is used,
%% e.g. if there are repeating distance entries, the order is like 1,2,3,3,4,..,n-1.
%%
%% @param dis is a symmetric distance matrix.
%%
%% @return disRank is the column rank matrices of X.
%%
function [disRank]=DistRanks(dis)

[n,m]=size(dis);
disRank=zeros(n,m);
for i=1:m
    [~,~,a]=unique(dis(:,i));
    disRank(:,i)=a;
end
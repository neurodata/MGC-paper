function [disRank]=disToRanks(dis)
% An auxiliary function that sorts the entries within each column by ascending order.
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
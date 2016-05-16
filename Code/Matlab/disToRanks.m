function [disRank]=disToRanks(dis)
% Transform from distance to ranking within each column, order from 1,...,n.
% For ties, the minimum ranking is used.

[n,m]=size(dis);
disRank=zeros(n,n);
for i=1:m
    [~,~,a]=unique(dis(:,i));
    disRank(:,i)=a;
%     % Alternative average rank
%     v=floor(tiedrank(dis(:,i)));
%     v(v==v(i))=1;
%     disRank(:,i)=v-1;
end
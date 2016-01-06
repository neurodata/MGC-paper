function [disRank]=disToRanks(dis)
% Transform from distance to ranking within each column, order from 0,...,n-1.
% For ties, the minimum ranking is used.

n=size(dis,1);
disRank=zeros(n,n);
for i=1:n
    [~,~,a]=unique(dis(:,i));
    disRank(:,i)=a-1;
    % Alternative average rank
    % v=floor(tiedrank(dis(:,i)));
    % v(v==v(i))=1;
    % disRank(:,i)=v-1;
end
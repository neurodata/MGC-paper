function [disRank]=disToRanks(dis)
% An auxiliary function of LocalCorr that calculates ranking within each column for a distance matrix, order from 1,...,n.
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
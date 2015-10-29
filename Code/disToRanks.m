function [disRank]=disToRanks(dis) %transform from distance to ranking, order from 0,...,n-1. 

n=size(dis,1);
[~,ind]=sort(dis(1:n,1:n));
disRank=zeros(n,n);
for i=1:n
    v=ind(1:n,i);
    for j=1:n
        disRank(v(j),i)=j-1;
    end
    if disRank(i,i)~=0 %change the diagonal rank to zero, if necessary
        disRank(v(1),i)=disRank(i,i);
        disRank(i,i)=0;
    end
end
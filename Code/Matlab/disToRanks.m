function [disRank]=disToRanks(dis)
% Transform from distance to ranking, order from 0,...,n-1.
% For ties, the minimum ranking is used.

n=size(dis,1);
disRank=zeros(n,n);
for i=1:n
    %                 v=floor(tiedrank(dis(:,i)));
    %                 v(v==v(i))=1;
    %                 disRank(:,i)=v-1;
    [~,~,a]=unique(dis(:,i));
    disRank(:,i)=a-1;
end

% n=size(dis,1);
% [~,ind]=sort(dis(1:n,1:n));
% disRank=zeros(n,n);
% for i=1:n
%     v=ind(1:n,i);
%     for j=1:n
%         disRank(v(j),i)=j-1;
%     end
%     if disRank(i,i)~=0 %change the diagonal rank to zero, if necessary
%         disRank(v(1),i)=disRank(i,i);
%         disRank(i,i)=0;
%     end
% end
function neighbor=verifyNeighbors(p,thres)
if nargin<2
    thres=0;
end
n=size(p,1);
pmin=min(min(p));

nind=find((p-pmin)<=thres);
neighbor=nind;
% neighbor=zeros(length(nind),1);
% 
% neighbor(:,2)=ceil(nind/n);
% neighbor(:,1)=nind-(neighbor(:,2)-1)*(n);

% ind=((p-pmin)<thres);
% nind=find(ind==1,1,'last');
% n=size(p,1);
% neighbor(2)=ceil(nind/n);
% neighbor(1)=nind-(neighbor(2)-1)*(n);

% if (abs(pmin-p(end,end))<thres)
%     pmin=p(end,end);
%     neighbor=size(p);
% else
%     neighbor(2)=find(min(p)==pmin,1,'last');
%     neighbor(1)=find(p(:,neighbor(2))==pmin,1,'last');
% end
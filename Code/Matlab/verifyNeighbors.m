function [neighbor]=verifyNeighbors(p,thres,all)
% return all indices whose value are no more than the minimal value by
% thres; if there are more than one such scale, return the largest possible
if nargin<2
    thres=0;
end
if nargin<3
    all=0; % if all is not 0, return all indices; if all is 0, return only the largest possible index.
end
n=size(p,1);
pmin=min(min(p));

nind=find((p-pmin)<=thres);
neighbor=nind(1);

if all==0
    max=0;
    for i=1:length(nind)
        l=ceil(nind(i)/n);
        k=nind(i)-(l-1)*n;
        if (k*l)>=max
            neighbor=nind(i);
            max=k*l;
        end
    end
end

    
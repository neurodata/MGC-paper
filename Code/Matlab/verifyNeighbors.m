function [neighbor]=verifyNeighbors(p,thres)
% return all indices whose value are no more than the minimal value by
% thres
if nargin<2
    thres=0;
end
n=size(p,1);
pmin=min(min(p));

nind=find((p-pmin)<=thres);
neighbor=nind;
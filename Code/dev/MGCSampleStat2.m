function [statMGC]=MGCSampleStat2(A,B,option)
% Author: Cencheng Shen
% This function estimate the Oracle MGC (i.e., optimal local correlation)
% from the local correlation map, which we call sample MGC statistic.
%
% It finds the largest connected region in the correlation map, such that
% each correlation is significant, i.e., larger than a certain threshold.
% To avoid correlation inflation by sample noise, it then computes Sample MGC as follows:
% for the largest correlation in the region, calculate the two minimal correlations
% along adjacent row scales and adjacent column scales, then take the larger one as the Sample MGC.
% If the region area is too small, or the estimated Sample MGC is no larger than the global, use the global correlation instead.
%
% Input: either a size m*n local correlation map, or two n*n distance matrices A and B and the global corr option.
% Output: the sample MGC statistic within [-1,1].
if nargin<2
    localCorr=A; % if there is only one input, asume the localCorr is given as A
else
    if nargin<3
        option='mcor';
    end
    localCorr=MGCLocalCorr(A,B,option); % otherwise compute the localCorr from given distance matrices
end
[m,n]=size(localCorr);
statMGC=localCorr(end); % take the global correlation by default
if m==1 || n==1;
    statMGC=localCorr(end);
    return;
end

% negCorr=localCorr(2:end,2:end);
% negCorr=negCorr(negCorr<0); % negative correlations
thres1=max(min(m,n),30);
thres2=2/thres1;
deg=(thres1*(thres1-3)/2-1)/2;
% tmp2=reshape((localCorr(2:m,2:n)+1)/2,1,(m-1)*(n-1));

p1=cdf('Beta',localCorr(end,end)/2+1/2,deg,deg);

% if p1>0.99
    thres1=statMGC+(icdf('Beta',0.95,deg,deg)*2-1);
% else
%     thres1=icdf('Beta',0.999,deg,deg)*2-1;
% end
% thres1=4*max(norm(negCorr,'fro')/sqrt(length(negCorr)),0.01); % threshold based on negative correlations
tau=0.1;

R=(localCorr>thres1); % find all correlations that are larger than the thresholds
localCorr(localCorr<=thres1)=0;
% R=(localCorr>0);
% % mean(mean(R))
% % find the largest connected component of all significant correlations
CC=bwconncomp(R,4);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
if isempty(idx)==false
    CC=CC.PixelIdxList{idx};
    R=zeros(m,n);
    R(CC)=1;
else
    R=0;
end
% if mean(mean(R))>=thres2
%     tmp=max(localCorr(R==1));
%     if tmp >= statMGC
%         statMGC=tmp;
%     end
% end
% thres2=0.2;
% if (mean(mean(R))>=thres2) && (max(max(localCorr))==max(localCorr(R==1)))
%     statMGC=max(statMGC,max(localCorr(R==1)));
% end
max1=max(max(localCorr));

% if (mean(mean(R))>=thres2) & (max1==max(localCorr(R==1)))
%     tmp=max(localCorr(R==1));
%     if tmp>statMGC;
%     statMGC=tmp;
%     end
% end
% % 
if (mean(mean(R))>thres2) && (max1==max(localCorr(R==1))) % proceed only when the region area is sufficiently large
    %[k,l]=find((localCorr>=max(localCorr(R==1)))&(R==1)); % find the scale within R that has the maximum correlation
    [k,l]=find(localCorr>=max1,1,'last');
    ln=ceil(tau*m); % boundary for checking adjacent rows
    km=ceil(tau*n); % boundary for checking adjacent columns
    for i=1:length(k)
        ki=k(i);
        li=l(i);
        
        % ensure the adjacent rows does not exceed the local correlation size, same for columns
        left=max(2,li-ln);
        right=min(n,li+ln);
        upper=max(2,ki-km);
        down=min(m,ki+km);
        
        tmp1=min(localCorr(upper:down,li)); % minimal correlation at given row and adjacent columns
        tmp2=min(localCorr(ki,left:right)); % minimal correlation at given column and adjacent rows
        tmp=max(tmp1,tmp2); % take the max of the two minimal correlations
        if tmp >= statMGC
            statMGC=tmp;
        end
    end
end

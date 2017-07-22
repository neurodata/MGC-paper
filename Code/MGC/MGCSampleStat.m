function [statMGC]=MGCSampleStat(A,B,option)
% Author: Cencheng Shen
% This function finds the sample MGC statistic by smoothing on the maximal local correlation, 
% which is an estimate of the Oracle MGC (i.e., optimal local correlation
% in terms of testing power).
%
% It first finds the largest connected region in the correlation map, such that
% each correlation is significant, i.e., larger than a certain threshold.
% To avoid correlation inflation by sample noise, it computes Sample MGC by smoothing the maximal correlation as follows:
% for the largest correlation within the region, calculate the two minimal correlations
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
statMGC=localCorr(end); % take the global correlation by default for sample mgc
if m==1 || n==1
    return;
end
prt=0.975;% percentile to consider as significant
tau=0.1; % percentage of adjacent scales to smooth with
  
% approximate a threshold for significant local dcorr
mn=max(min([m,n,80]),20); 
% The degree mn equals sample size, but is otherwise bounded by [20,80]: less than 20 is too small for the
% approximation to be accurate, while larger than 80 caused the threshold to be too small
if mn>20
    % for sufficient sample size, estimate based on normal distribution approximation from Szekely2013
    thres=sqrt(mn*(mn-3)/2-1);
    thres=icdf('normal',prt,0,1)/thres;
    %%% thres=mn*(mn-3)/4-1/2;
    %%% thres=icdf('beta',prtl,thres,thres)*2-1;
else
    % for insufficient sample size, estimate based on nonparamtric estimation via negative statistics
    thres=localCorr(2:end,2:end);
    thres=thres(thres<0); % negative correlations
    thres=3.5*max(norm(thres,'fro')/sqrt(length(thres)),0.01);
end

% find all correlations that are larger than the threshold
%localCorr(localCorr<=thres)=0; % localCorr(localCorr<=max(thres,thres2))=0;
R=(localCorr>thres);
% find the largest connected component of all significant correlations
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

% Smooth maximal local correlation in the significant region for sample mgc
if mean(mean(R))>=1/mn % proceed only when the significant region is sufficiently large
    [k,l]=find((localCorr>=max(localCorr(R==1)))&(R==1)); % find all scales within R that maximize the local correlation
    
    ln=ceil(tau*m); % number of adjacent rows to check
    km=ceil(tau*n); % number of adjacent columns to check
    for i=1:length(k)
        ki=k(i);
        li=l(i);
        
        % boundary of rows and columns for smoothing
        left=max(2,li-ln);
        right=min(n,li+ln);
        upper=max(2,ki-km);
        down=min(m,ki+km);
        
        tmp1=min(localCorr(upper:down,li)); % take minimal correlation at given row and along adjacent columns
        tmp2=min(localCorr(ki,left:right)); % take minimal correlation at given column and along adjacent rows
        tmp=max(tmp1,tmp2); % take the max of the two minimal correlations for sample mgc
        
        % use the smoothed maximal local correlation for mgc statistic, only when it is larger than global dcorr
        if tmp >= statMGC 
            statMGC=tmp;
        end
    end
end

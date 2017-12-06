%% The main function that computes the MGC measure between two datasets:
%% It first computes all local correlations,
%% then use the maximal statistic among all local correlations based on thresholding.
%%
%% @param A is a distance matrix or a n*d data matrix; if it is not a square matrix with zeros on diagonal, it is treated as n*d data;
%% @param B is a second distance matrix or a n*d data matrix, with the same distance matrix check as A;
%% @param option is a string that specifies which global correlation to build up-on, including 'mgc','dcor','mantel', and 'rank'.
%%
%% @return A list contains the following output:
%% @return statMGC is the sample MGC statistic within [-1,1];
%% @return localCorr consists of all local correlations by double matrix index;
%% @return optimalScale the estimated optimal scale in matrix single index.
%%
%% @export
%%
function [statMGC, localCorr, optimalScale]=MGCSampleStat(A,B,option)
if nargin<3
    option='mgc';
end
localCorr=MGCLocalCorr(A,B,option); % compute all localCorr
[m,n]=size(localCorr);
if m==1 || n==1
    statMGC=localCorr(end);
    optimalScale=m*n;
else
    sz=size(A,1)-1; % sample size minus one
    R=Thresholding(localCorr,m,n,sz); % find a connected region of significant local correlations
    [statMGC,optimalScale]=Smoothing(localCorr,m,n,R); % find the maximal within the significant region
end

%% An auxiliary function that finds a region of significance in the local correlation map by thresholding.
%%
%% @param localCorr is all local correlations;
%% @param m is the number of rows of localCorr;
%% @param n is the number of columns of localCorr;
%% @param sz is the sample size of original data (which may not equal m or n in case of repeating data).
%%
%% @return R is a binary matrix of size m and n, with 1's indicating the significant region.
%%
function R=Thresholding(localCorr,m,n,sz)
% A threshold is estimated based on normal distribution approximation from Szekely2013
prt=1-0.02/sz; % percentile to consider as significant
%thres=sqrt(sz*(sz-3)/2-1); % normal approximation, which is equivalent to beta approximation for n larger than 10
%thres=icdf('normal',prt,0,1)/thres;
thres=sz*(sz-3)/4-1/2; % beta approximation
thres=(betainv(prt,thres,thres))*2-1;

opt=0; % set opt=1 and add the negative local correlation as a non-parametric and data-adaptive threshold
if opt==1
    thres1=localCorr;
    thres1=thres1(thres1<0); % all negative correlations
    thres1=5*norm(thres1,'fro')/sqrt(length(thres1)); % 5 times the standard deviation of negative correlations
    thres=max(thres1,thres); % Use the maximal of paratemetric and non-parametric thresholds
end
thres=max(thres,localCorr(end)); % take the maximal of threshold and local correlation at the maximal scale

% Find the largest connected component of significant correlations
R=(localCorr>thres);
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
% Find all correlations that are larger than the threshold
% [~,~,~,R]=FindLargestRectangles((localCorr>=thres), [0 0 1],[2,2]);
% % optimalInd=find(optimalInd==1);
% R=double(R);
% if sum(sum(R))==0
%     R=0;
% end

%% An auxiliary function that finds the maximal within the significant region R:
%% If area of R is too small, return the last local corr; otherwise take the maximum within R.
%%
%% @param localCorr is all local correlations;
%% @param m is the number of rows of localCorr;
%% @param n is the number of columns of localCorr;
%% @param R is a binary matrix of size m by n indicating the significant region.
%%
%% @return A list contains the following output:
%% @return statMGC is the sample MGC statistic within [-1,1];
%% @return optimalScale the estimated optimal scale in matrix single index.
%%
function [statMGC, optimalScale]=Smoothing(localCorr,m,n,R)
statMGC=localCorr(end); % default sample mgc to local corr at maximal scale
optimalScale=m*n; % default the optimal scale to maximal scale
if (norm(R,'fro')~=0)
    % tau=0; % number of adjacent scales to smooth with
    if sum(sum(R))>=2*min(m,n) % proceed only when the region area is sufficiently large
        tmp=max(localCorr(R==1));
        [k,l]=find((localCorr>=tmp)&(R==1)); % find all scales within R that maximize the local correlation
        
%         ln=ceil(tau); % number of adjacent rows to check
%         km=ceil(tau); % number of adjacent columns to check
%         for i=1:length(k)
%             ki=k(i);
%             li=l(i);
%             
%             % index of adjacent rows and columns
%             left=max(2,li-ln);
%             right=min(n,li+ln);
%             upper=max(2,ki-km);
%             down=min(m,ki+km);
%             
%             tmp1=min(localCorr(upper:down,li)); % take minimal correlation at given row and along adjacent columns
%             tmp2=min(localCorr(ki,left:right)); % take minimal correlation at given column and along adjacent rows
%             tmp=max(tmp1,tmp2); % take the min for sample mgc
            
            if tmp >= statMGC
                statMGC=tmp;
                optimalScale=(l-1)*m+k; % take the scale of maximal stat and change to single index
            end
%         end
    end
end
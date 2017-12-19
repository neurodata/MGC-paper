%% The main function that tests independent between two data sets by MGC and permutation test.
%%
%% @param X is a distance matrix or a n*d data matrix; if it is not a square matrix with zeros on diagonal, it is treated as n*d data; 
%% @param Y is a second distance matrix or a n*d data matrix, with the same distance matrix check as X;
%% @param rep specifies the number of replicates to use for the permutation test;
%% @param option is a string that specifies which global correlation to build up-on, including 'mgc','dcor','mantel', and 'rank'.
%%
%% @return A list contains the following output:
%% @return P-value and test statistic of MGC;
%% @return All local p-values by double matrix index, all local correlations by double matrix index, and the estimated optimal scales as matrix single indices.
%%
%% Note that one should avoid report positive discovery via minimizing individual p-values of local correlations,
%% unless corrected for multiple testing problem.
%%
%% @export
%% 
function  [pMGC,statMGC,pLocalCorr,localCorr,optimalScale]=MGCPermutationTest(X,Y,rep,option)

if nargin<3
    rep=1000; % use 1000 random permutations by default
end
if nargin<4
    option='mgc';  % use mgc by default
end
% Use the data size and diagonal element to determine if the given data is a distance matrix or not
if size(X,1)~=size(X,2) || sum(diag(X).^2)>0
    X=squareform(pdist(X));
%     disp('The first data is not a Euclidean distance matrix; transformed to distance matrix instead.')
end
if size(Y,1)~=size(Y,2) || sum(diag(Y).^2)>0
    Y=squareform(pdist(Y));
%     disp('The second data is not a Euclidean distance matrix; transformed to distance matrix instead.')
end
np=size(Y,1);

% Compute sample MGC, all local correlations, and the estimated optimal scale
[statMGC,localCorr, optimalScale]=MGCSampleStat(X,Y,option);
[m,n]=size(localCorr);
pLocalCorr=zeros(m,n);pMGC=0;

% Compute sample MGC and all local correlations for each permuted data
for r=1:rep
    % Use random permutations on the second data set
    per=randperm(np);
    YN=Y(per,per);
    [tmp2,tmp]=MGCSampleStat(X,YN,option);
    pMGC=pMGC+(tmp2>=statMGC)/rep;
    pLocalCorr=pLocalCorr+(tmp>=localCorr)/rep;
end

% % Estimate the optimal scales via both the statistics and p-values if possible
% [~,~,~,optimalScale]=FindLargestRectangles((pLocalCorr<=pMGC)&(localCorr>=statMGC), [0 0 1]);
% optimalScale=find(optimalScale==1);
% if (isempty(optimalScale))
%     [~,~,~,optimalScale]=FindLargestRectangles((pLocalCorr<=pMGC), [0 0 1]);
%     optimalScale=find(optimalScale==1);
%     if (isempty(optimalScale))
%         optimalScale=m*n; % if the global scale is not selected in the largest rectangle while being optimal, we take the global scale instead.
%     end
% end
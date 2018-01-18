function  [p, test]=PermutationTest2(C,D,rep,option)
% Author: Cencheng Shen
% This is for HHG permutation test.
% The outputs are the p-values of HHG and the test statistic.
if nargin<3
    rep=1000;
end
n=size(C,1);

% Calculate the observed test statistics for the given data sets
p=0;
if strcmp(option,'corr')
    ch=1;
else
    if strcmp(option,'HSIC')
        ch=2;     
    else
        ch=3;   
    end
end

switch ch
    case 1
        test=corr(C,D);
    case 2
        test=HSIC(C,D);
    case 3
        an=mine(C',D');
        test=an.mic;
end
% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % Use random permutations;
    per=randperm(n);
    DN=D(per,:);
    switch ch
        case 1
            tmp=corr(C,DN);
        case 2
            tmp=HSIC(C,DN);
        case 3
            an=mine(C',DN');
            tmp=an.mic;
    end
    p=p+(tmp>=test)/rep;
end
% if (sum(pHHG<1/rep)>0)
%     pHHG=pHHG+1/rep;
% end
function [p1,p2,p3, p4,p5,p6,p7,p1All,p2All,p3All]=CorrPermDistTest(C,D,rep,titlechar,option)
% Author: Cencheng Shen
% Permutation Tests for identifying dependency.
% The output are the p-values of MGC by dcorr/mcorr/Mantel, and global dcorr/mcorr/Mantel/HHG.

% Parameters:
% C & D should be two n*n distance matrices,
% rep specifies the number of random permutations to use for the permutation test,
% alpha specifies the type 1 error level,
% option specifies whether each test statistic is calculated or not,
% neighborhood can be specified beforehand so as to skip the optimal scale
% estimation.
if nargin<4
    titlechar='RealData';
end
if nargin<5
    option=[1,2,3,4]; % Default option. Setting any to 0 to disable the calculation of MGC by dcorr/mcorr/Mantel, global dcorr/mcorr/Mantel, HHG, in order; set the first three
end

% Global Correlations
p1All=1;p2All=1;p3All=1;p1=1;p2=1;p3=1;ind1=0;ind2=0;ind3=0;p7=0;t1All=0;t2All=0;t3All=0;
if option(1)==1
    [p1,p1All,~,t1All,ind1]=MGCPermutationTest(C,D,rep,option(1));
end
if option(2)==2
    [p2,p2All,~,t2All,ind2]=MGCPermutationTest(C,D,rep,option(2));
end
if option(3)==3
    [p3,p3All,~,t3All,ind3]=MGCPermutationTest(C,D,rep,option(3));
end
if option(4)==4
    [p7]=HHGPermutationTest(C,D,rep);
end
p4=p1All(end);p5=p2All(end);p6=p3All(end);

% if p2<0.05 || p5<0.05
%     p2
%     p5
% end
%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-3));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

pre1=strcat(rootDir,'Data/Results/');% The folder to save figures
filename=strcat(pre1,'CorrPermDistTestType',titlechar);
save(filename,'titlechar','rep','option','p1All','p2All','p3All','p4','p5','p6','p7','p1','p2','p3','ind1','ind2','ind3','t1All','t2All','t3All');

function  [p, test]=HHGPermutationTest(C,D,rep)
% Author: Cencheng Shen
% This is for HHG permutation test.
% The outputs are the p-values of HHG and the test statistic.
if nargin<3
    rep=1000;
end
n=size(C,1);

% Calculate the observed test statistics for the given data sets
p=0;
test=HHG(C,D);

% Now Permute the second dataset for rep times, and calculate the p-values
for r=1:rep
    % Use random permutations;
    per=randperm(n);
    DN=D(per,per);
    tmp=HHG(C,DN);
    p=p+(tmp<test)/rep;
end
% Output the p-value
p=1-p;
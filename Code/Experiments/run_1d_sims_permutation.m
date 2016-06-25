function p1=run_1d_sims_permutation(rep1,rep2)
% run 1-dimensional simulations based on MGC approximate scale for permutation test

%%%
fpath = mfilename('fullpath');
findex=strfind(fpath,'/');
if isempty(findex)
    findex=strfind(fpath,'\');
end
rootDir=fpath(1:findex(end-2));
p = genpath(rootDir);
gits=strfind(p,'.git');
colons=strfind(p,':');
for i=0:length(gits)-1
    endGit=find(colons>gits(end-i),1);
    p(colons(endGit-1):colons(endGit)-1)=[];
end
addpath(p);

if nargin < 1
    rep1=1000; % number of MC replicates for power computation
end
if nargin < 2
    rep2=1000; % number of random permutations for p-value calculation
end
%%
thres=0.8;dim=1;noise=1;type=1:20;
[p1]=CorrSimPermScale1(type,dim,thres,rep1,rep2,noise);
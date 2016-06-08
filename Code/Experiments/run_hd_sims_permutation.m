function p1=run_hd_sims_permutation(rep1,rep2)
% run high-dimensional simulations based on MGC approximate scale for permutation test
if nargin < 1
    rep1=1000; % number of MC replicates for power computation
end
if nargin < 2
    rep2=1000; % number of random permutations for p-value calculation
end
%%
thres=0.5;dim=2;noise=0;type=1:20;
[p1]=CorrSimPermScale1(type,dim,thres,rep1,rep2,noise);
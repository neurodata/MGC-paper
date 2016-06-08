function run_outlier_sims(rep1,rep2)
% run the outlier model simulations
if nargin < 1
    rep1=200; % number of MC replicates for MGC scale estimation
end
if nargin < 2
    rep2=500; % number of MC replicates for power computation
end

n=100; dim=1; lim=1; 
noise=0.3;
CorrIndTest(0,n,dim,lim,rep1,rep2,noise);
noise=0.5;
CorrIndTest(0,n,dim,lim,rep1,rep2,noise);
noise=0.7;
CorrIndTest(0,n,dim,lim,rep1,rep2,noise);
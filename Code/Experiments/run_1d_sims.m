function run_1d_sims(rep,type)
% run 1-dimensional simulations
% run_1d_sims(2000,10000) % the replicates for draft

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));

if nargin < 1
    rep=100; % number of MC replicates for MGC scale estimation
end
if nargin < 2
    type=1:20; 
end

%Ind
n=100;lim=20;dim=1;
for i=type
    CorrIndTest(i,n,dim,lim,rep); % the output are saved to ../../data/results
end
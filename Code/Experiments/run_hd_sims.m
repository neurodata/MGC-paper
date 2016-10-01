function run_hd_sims(rep1,rep2,type)
% run high-dimensional simulations
% the output are saved to ../../data/results

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

if nargin < 1
    rep1=2000; % number of MC replicates for MGC scale estimation
end
if nargin < 2
    rep2=10000; % number of MC replicates for power computation
end
if nargin < 3
    type=1:20;
end

n=100;
for i=type
    switch i
        case {1,2,3}
            dim=1000;
        case {5,6,7,15,11,16,17,8}
            dim=20;
        case {4,12,13,19}
            dim=10;
        case {14,18}
            dim=40;
        case {10,9,20}
            dim=100;
    end
    if dim<20
        lim=10;
    else
        lim=20;
    end
    CorrIndTestDim(i,n,dim,lim,rep1,rep2);
end
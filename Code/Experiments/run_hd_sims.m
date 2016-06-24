function run_hd_sims(rep1,rep2)
% run high-dimensional simulations
% the output are saved to ../../data/results

%%
fpath = mfilename('fullpath');
findex=strfind(fpath,'\');
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
    rep1=2000; % number of MC replicates for MGC scale estimation
end
if nargin < 2
    rep2=10000; % number of MC replicates for power computation
end

n=100; dim=1000; lim=20; 
CorrIndTestDim(1,n,dim,lim,rep1,rep2);
CorrIndTestDim(2,n,dim,lim,rep1,rep2);
CorrIndTestDim(3,n,dim,lim,rep1,rep2);
dim=20;lim=20;
CorrIndTestDim(4,n,dim,lim,rep1,rep2);
dim=10;lim=10;
CorrIndTestDim(5,n,dim,lim,rep1,rep2);
dim=40;lim=20;
CorrIndTestDim(6,n,dim,lim,rep1,rep2);
CorrIndTestDim(7,n,dim,lim,rep1,rep2);
CorrIndTestDim(8,n,dim,lim,rep1,rep2);
dim=40;lim=20;
CorrIndTestDim(9,n,dim,lim,rep1,rep2);
dim=100;
CorrIndTestDim(10,n,dim,lim,rep1,rep2);


n=100; dim=20; lim=20; 
CorrIndTestDim(11,n,dim,lim,rep1,rep2);
CorrIndTestDim(12,n,dim,lim,rep1,rep2);
CorrIndTestDim(13,n,dim,lim,rep1,rep2);
dim=40; lim=20;
CorrIndTestDim(14,n,dim,lim,rep1,rep2);
CorrIndTestDim(15,n,dim,lim,rep1,rep2);
dim=10;lim=10;
CorrIndTestDim(16,n,dim,lim,rep1,rep2);
dim=10;lim=10;
CorrIndTestDim(17,n,dim,lim,rep1,rep2);
dim=1000;lim=20;
CorrIndTestDim(18,n,dim,lim,rep1,rep2);
dim=100;
CorrIndTestDim(19,n,dim,lim,rep1,rep2);
CorrIndTestDim(20,n,dim,lim,rep1,rep2);
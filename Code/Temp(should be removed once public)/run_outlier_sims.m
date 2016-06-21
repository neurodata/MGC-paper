function [p1,p2]=run_outlier_sims(rep1,rep2)
% run the outlier model simulations
if nargin < 1
    rep1=200; % number of MC replicates for MGC scale estimation
end
if nargin < 2
    rep2=500; % number of MC replicates for power computation
end

% n=100; dim=1; lim=1; 
% noise=0.3;
% CorrIndTest(0,n,dim,lim,rep1,rep2,noise);
% noise=0.5;
% CorrIndTest(0,n,dim,lim,rep1,rep2,noise);
% noise=0.7;
% CorrIndTest(0,n,dim,lim,rep1,rep2,noise);

n=100;dim=1;rep1=100;rep2=100;noise=0.1;
p1=zeros(n,n);
p2=0;
alpha=0.05;
for r1=1:rep1
        [x, y]=CorrSampleGenerator(0,n,dim,1, noise);
        C=squareform(pdist(x));
        D=squareform(pdist(y));
        [~,~,~,~,~,~,~,~,pp2]=CorrPermDistTest(C,D,rep2,'PermInd',[0;2;0;0]);
        p1=p1+(pp2<alpha)/rep1;
        p2=p2+(pp2(end)<alpha)/rep1;
%         p
    end
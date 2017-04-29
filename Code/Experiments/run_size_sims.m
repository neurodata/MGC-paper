function run_size_sims(rep,dim,type)
% run 1-dimensional simulations
% run_size_sims(100,10,1:19)
% run_size_sims(100,1,1:19)
%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));

pre1=strcat(rootDir,'Data/Results/');% The folder to save figures
filename=strcat(pre1,'CorrIndSize','Dim',num2str(dim));

if nargin < 1
    rep=200; % number of MC replicates for MGC scale estimation
end
if nargin < 2
    dim=1;
end
if nargin < 3
    type=1:19; 
end

nMax=1000;
thres=0.85;
if dim>1
    noise=0;
else
    noise=1;
end
alpha=0.05;
%Ind
SampleSize=zeros(6,max(type)+1);
for i=type
    for option=0:5
        SampleSize(option+1,i)=CorrIndTestSize(i,nMax,dim,thres,rep,noise,alpha,option);
        if option==0 && SampleSize(1,i)==1000
            SampleSize(2:end,i)=1000;
            break;
        end
    end
    save(filename,'SampleSize','nMax','dim','thres','rep','noise','alpha','type');
end

SampleSize(:,max(type)+1)=mean(SampleSize(:,type),2);
SampleSize=floor(SampleSize);

% Save the results
%%% File path searching
save(filename,'SampleSize','nMax','dim','thres','rep','noise','alpha','type');

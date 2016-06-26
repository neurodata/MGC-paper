function run_fig1Data(type,n,noise)
% generate data for figure 1
%%% File path searching
if nargin<1
    type=11;
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

%type=11;
option=2;
%n=50;
dim=1;
%noise=0;
rep=1000;

% optimal scale
tA=zeros(n,n,rep);
tN=zeros(n,n,rep);
    for r=1:rep;
        [x, y]=CorrSampleGenerator(type,n,dim,1, noise);
        CA=squareform(pdist(x));
        DA=squareform(pdist(y));
        tA(:,:,r)=LocalCorr(CA,DA,2);
        [x, y]=CorrSampleGenerator(type,n,dim,0, noise);
        CA=squareform(pdist(x));
        DA=squareform(pdist(y));
        tN(:,:,r)=LocalCorr(CA,DA,2);
    end
    power1=zeros(n,n);
    alpha=0.05;
    for i=1:n;
        for j=1:n;
            dCorT=sort(tN(i,j,:),'descend');
            cut1=dCorT(ceil(rep*alpha));
            power1(i,j)=mean(tA(i,j,:)>cut1);
        end
    end
power1(1,:)=0;power1(:,1)=0; % Set the powers of all local tests at rank 0 to 0
neighbor=maxNeighbors(power1,0,tN,tA);


% Generate data
[x, y]=CorrSampleGenerator(type,n,dim,1, noise);
if noise~=0
    [x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
end
% Get distance matrix
[x,ind]=sort(x,'ascend');
y=y(ind);
C=squareform(pdist(x));
D=squareform(pdist(y));

% Permutation p-value
tA=LocalCorr(C,D,2);
tN=zeros(rep,n,n);
pAll=zeros(n,n);
for r=1:rep;
    per=randperm(n);
    tmp=LocalCorr(C,D(per,per),2);
    tN(r,:,:)=tmp;
    if r==1
        pAll=(tmp<tA)/rep;
    else
        pAll=pAll+(tmp<tA)/rep;
    end
end

l=ceil(neighbor/n);
k=neighbor-n*(l-1);

save(strcat(rootDir,'Data/Results/CorrFigure1Type',num2str(type),'.mat'),'tA','tN','type','n','option','dim','noise','rep','power1','neighbor','pAll','k','l','C','D','x','y');
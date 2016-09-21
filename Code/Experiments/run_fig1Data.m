function run_fig1Data(type,n,noise,rep)
% generate data for figure 1
%%% File path searching
if nargin<1
    type=8;
end
if nargin<2
    n=50;
end
if nargin<3
    noise=0;
end
if nargin<4
    rep=1000;
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));

%type=11;
%n=50;
dim=1;
%noise=0;

% optimal scale
tA=zeros(n,n,rep);
tN=zeros(n,n,rep);
for r=1:rep;
    [x, y]=CorrSampleGenerator(type,n,dim,1, noise);
    CA=squareform(pdist(x));
    DA=squareform(pdist(y));
    tA(:,:,r)=MGCLocalCorr(CA,DA,'mcor');
    [x, y]=CorrSampleGenerator(type,n,dim,0, noise);
    CA=squareform(pdist(x));
    DA=squareform(pdist(y));
    tN(:,:,r)=MGCLocalCorr(CA,DA,'mcor');
end
powerMLocal=zeros(n,n);
alpha=0.05;
for i=1:n;
    for j=1:n;
        dCorT=sort(tN(i,j,:),'descend');
        cut1=dCorT(ceil(rep*alpha));
        powerMLocal(i,j)=mean(tA(i,j,:)>cut1);
    end
end
powerMLocal(1,:)=0;powerMLocal(:,1)=0; % Set the powers of all local tests at rank 0 to 0
neighbor=maxNeighbors(powerMLocal,0,tN,tA);


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
tA=MGCLocalCorr(C,D,'mcor');
test=MGCSampleStat(tA);
tN=zeros(rep,n,n);
testN=zeros(rep,1);
pMLocal=zeros(n,n);
pMGC=0;
for r=1:rep;
    per=randperm(n);
    tmp=MGCLocalCorr(C,D(per,per),'mcor');
    tmp2=MGCSampleStat(tmp);
    tN(r,:,:)=tmp;
    testN(r)=tmp2;
    pMLocal=pMLocal+(tmp>=tA)/rep;
    pMGC=pMGC+(tmp2>=test)/rep;
end
% set the p-values of local corr at rank 0 to maximum, since they should not be used
if min(min(pMLocal(2:end,2:end)))==0
    pMLocal=pMLocal+1/rep;
end
if min(min(pMLocal(2:end,2:end)))>pMGC
    pMGC=min(min(pMLocal(2:end,2:end)));
end
pMLocal(pMLocal>1)=1;
pMLocal(1,:)=1;pMLocal(:,1)=1;

% find optimal scale
warning('off','all');
[~,~,~,optimalInd]=FindLargestRectangles((pMLocal<=pMGC), [0 0 1],[2,2]);
optimalInd=find(optimalInd==1);
if (pMLocal(end)<pMGC && isempty(find(optimalInd==n*n, 1)))
    optimalInd=n*n;
end
l=ceil(neighbor/n);
k=neighbor-n*(l-1);

% Mcorr
RC=DistRanks(C);
RD=DistRanks(D)';

H=eye(n)-(1/n)*ones(n,n);
A=H*C-C/n;
B=D*H-D/n;
for j=1:n
    A(j,j)=0;
    B(j,j)=0;
end
mcorrH=A.*B;

% Mantel
C=C-sum(sum(C))/n/(n-1);
D=D-sum(sum(D))/n/(n-1);
mantelH=C.*D;

% MGC
RC=(RC>k);
RD=(RD>l);
A_MGC=A;
B_MGC=B;
A_MGC(RC)=0;
B_MGC(RD)=0;
% A_MGC=A_MGC;
% B_MGC=B_MGC;
%A_MGC(1:n+1:n^2)=0;
%B_MGC(1:n+1:n^2)=0;
C_MGC=(A_MGC-mean(mean(A_MGC))).*(B_MGC-mean(mean(B_MGC)));

save(strcat(rootDir,'Data/Results/CorrFigure1Type',num2str(type),'n',num2str(n),'.mat'),'tA','tN','test','testN','type','n','dim','noise','rep','powerMLocal','neighbor','pMLocal','pMGC','optimalInd','k','l','C','D','x','y','A','B','mantelH','mcorrH','A_MGC','B_MGC','C_MGC');
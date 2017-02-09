function run_realData(rep,select)
% run real data experiments
if nargin<1
    rep=10000;
end
if nargin<2
    select=4;
end

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));

%%%
if select==1 || select==4
    load(strcat(rootDir,'Data/Preprocessed/BrainCP.mat'))
    n=42;
    CorrPermDistTest(distC,distP,rep,'BrainCxP');
end

%%%
if select==2 || select==4
    load(strcat(rootDir,'Data/Preprocessed/BrainHippoShape.mat'))
    n=114;
    y=squareform(pdist((Label+unifrnd(0,0.01,n,1))));
    %     y=squareform(pdist(Label));
    %     y=y+1;
    %     for i=1:n
    %         y(i,i)=0;
    %     end
    %     CorrPermDistTest(LMLS,y,rep,'BrainLMLxY');
    CorrPermDistTest(LMRS,y,rep, 'BrainLMRxY');
    %     CorrPermDistTest(LMLS,LMRS,rep,'BrainLMLxLMR');
end

%%%
if select==3 || select==4
    load(strcat(rootDir,'Data/Preprocessed/semipar.mat'))
    n=109;
    distCCI=squareform(pdist(cci));
    CorrPermDistTest(distMigrain(ind,ind),distCCI(ind,ind),rep,'MigrainxCCI');
    % CorrPermDistTest(distM2g(ind,ind),distCCI(ind,ind),rep,'M2gxCCI');
    % CorrPermDistTest(distM2g(ind,ind),distMigrain(ind,ind),rep,'M2gxMigrain');
end
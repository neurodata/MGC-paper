function run_realData(rep,select)
% run real data experiments
if nargin<1
    rep=10000;
end
if nargin<2
    select=5;
end

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));

%%%
if select==1 || select==5
    load(strcat(rootDir,'Data/Preprocessed/BrainCP.mat'))
    n=42;
    CorrPermDistTest(distC,distP,rep,'BrainCxP');
end

%%%
if select==2 || select==5
    load(strcat(rootDir,'Data/Preprocessed/BrainHippoShape.mat'))
    n=114;
    y=squareform(pdist((Label+unifrnd(0,0.01,n,1))));
        %y=squareform(pdist(Label));
    %     y=y+1;
    %     for i=1:n
    %         y(i,i)=0;
    %     end
    %     CorrPermDistTest(LMLS,y,rep,'BrainLMLxY');
    CorrPermDistTest(LMRS,y,rep, 'BrainLMRxY');
    %     CorrPermDistTest(LMLS,LMRS,rep,'BrainLMLxLMR');
end

%%%
if select==3 || select==5
    load(strcat(rootDir,'Data/Preprocessed/semipar.mat'))
    n=109;
    distCCI=squareform(pdist(cci));
    CorrPermDistTest(distMigrain(ind,ind),distCCI(ind,ind),rep,'MigrainxCCI');
    % CorrPermDistTest(distM2g(ind,ind),distCCI(ind,ind),rep,'M2gxCCI');
    % CorrPermDistTest(distM2g(ind,ind),distMigrain(ind,ind),rep,'M2gxMigrain');
end

if select==4 || select==5
    load(strcat(rootDir,'Data/Preprocessed/proteomics.mat'))
    C=squareform(pdist(A'));
    D=squareform(pdist(LabelInd));
%     per1=find(LabelInd==1);
%     per2=find(LabelInd==2);
%     per3=find(LabelInd==3);
%     per4=find(LabelInd==4);
%     per5=find(LabelInd==5);
%     n1=size(per1,1);
%     n2=size(per2,1);
%     n3=size(per3,1);
%     n4=size(per4,1);
%     n5=size(per5,1);
%     t1=randperm(n1);
%     t2=randperm(n2);
%     t3=randperm(n3);
%     t4=randperm(n4);
%     t5=randperm(n5);
%     sub=5;
%     per=[per1(t1(1:sub));per2(t2(1:sub));per3(t3(1:sub));per4(t4(1:sub));per5(t5(1:sub))];
    per=(LabelInd>2);
    [pMGC,pD,pM,pP, pHHG]=CorrPermDistTest(C,D,rep,'ProtecomicsFull');
%     C1=squareform(pdist(A(:,per1)));
%     C2=squareform(pdist(A(:,per2)));
%     C3=squareform(pdist(A(:,per3)));
%     C4=squareform(pdist(A(:,per4)));
%     C5=squareform(pdist(A(:,per5)));
    [pMGC,pD,pM,pP, pHHG]=CorrPermDistTest(C(per,per),D(per,per),rep,'ProtecomicsPartial');
end
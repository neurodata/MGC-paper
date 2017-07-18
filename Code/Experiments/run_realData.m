function run_realData(rep,select)
% run real data experiments
% run_realData(10000)  % the replicates for draft
if nargin<1
    rep=100;
end
if nargin<2
    select=7;
end

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));

%%%
if select==1 || select==7
    load(strcat(rootDir,'Data/Preprocessed/BrainCP.mat'))
    n=42;
    CorrPermDistTest(distC,distP,rep,'BrainCxP');
end

%%%
if select==2 || select==7
    load(strcat(rootDir,'Data/Preprocessed/BrainHippoShape.mat'))
    n=114;
    y=squareform(pdist((Label+unifrnd(0,0.01,n,1))));
    %y=squareform(pdist(Label));
    %y=y+1;
    %for i=1:n
    %    y(i,i)=0;
    %end
    %CorrPermDistTest(LMLS,y,rep,'BrainLMLxY');
    CorrPermDistTest(LMRS,y,rep, 'BrainLMRxY');
    %CorrPermDistTest(LMLS,LMRS,rep,'BrainLMLxLMR');
end

%%%
if select==3 || select==7
    load(strcat(rootDir,'Data/Preprocessed/semipar.mat'))
    n=109;
    distCCI=squareform(pdist(cci));
    CorrPermDistTest(distMigrain(ind,ind),distCCI(ind,ind),rep,'MigrainxCCI');
    % CorrPermDistTest(distM2g(ind,ind),distCCI(ind,ind),rep,'M2gxCCI');
    % CorrPermDistTest(distM2g(ind,ind),distMigrain(ind,ind),rep,'M2gxMigrain');
end

if select==4 || select==7
    load(strcat(rootDir,'Data/Preprocessed/proteomics.mat'))
    
    % Ovarian vs Normal
    per=(LabelIndAll==1 | LabelIndAll==4);
    D=LabelIndAll(per);
    D=squareform(pdist(D));
    D(D>0)=1;
    m=318;
    C=squareform(pdist(A(:,per)'));
    CorrPermDistTest(C,D,rep,'ProteomicsOvavsNormal');
    
    % Screening Ovariance vs Normal
    pMGC=zeros(m,1);pD=zeros(m,1);pM=zeros(m,1);pP=zeros(m,1);pHHG=zeros(m,1);testMGC=zeros(m,1);testD=zeros(m,1);testM=zeros(m,1);testHHG=zeros(m,1);
    for i=1:m
        i
        C=squareform(pdist(A(i,per)'));
        [pMGC(i),pD(i),pM(i),pP(i), pHHG(i),testMGC(i),testD(i),testM(i),testHHG(i)]=CorrPermDistTest(C,D,rep);
    end
    
    testMGC(:,6)=pMGC;
    testMGC(:,7)=pP;
    testMGC(:,8)=pD;
    testMGC(:,9)=pM;
    testMGC(:,10)=pHHG;
    testMGC(:,2)=testD;
    testMGC(:,3)=testD;
    testMGC(:,4)=testM;
    testMGC(:,5)=testHHG;
    save(strcat(rootDir,'Data/Results/','ScreeningOvavsNormal'),'testMGC');
end

if select==5 || select==7
    load(strcat(rootDir,'Data/Preprocessed/proteomics.mat'))
    per=(LabelIndAll==1 | LabelIndAll==2);
    D=LabelIndAll(per);
    D=squareform(pdist(D));
    D(D>0)=1;
    m=318;
    
    C=squareform(pdist(A(:,per)'));
    CorrPermDistTest(C,D,rep,'ProteomicsPancvsNormal');
    % ind=181;%296 for hhg
    % C=squareform(pdist(A(ind,per)'));
    % CorrPermDistTest(C,D,rep*2,'ProteomicsPancvsNormalNeuroganin');
    
    % Screening Panc vs Normal
    pMGC=zeros(m,1);pD=zeros(m,1);pM=zeros(m,1);pP=zeros(m,1);pHHG=zeros(m,1);testMGC=zeros(m,1);testD=zeros(m,1);testM=zeros(m,1);testHHG=zeros(m,1);
    for i=1:m
        C=squareform(pdist(A(i,per)'));
        [pMGC(i),pD(i),pM(i),pP(i), pHHG(i),testMGC(i),testD(i),testM(i),testHHG(i)]=CorrPermDistTest(C,D,rep);
    end
    
    testMGC(:,6)=pMGC;
    testMGC(:,7)=pP;
    testMGC(:,8)=pD;
    testMGC(:,9)=pM;
    testMGC(:,10)=pHHG;
    testMGC(:,2)=testD;
    testMGC(:,3)=testD;
    testMGC(:,4)=testM;
    testMGC(:,5)=testHHG;
    save(strcat(rootDir,'Data/Results/','ScreeningPancvsNormal'),'testMGC');
end


if select==6 || select==7
    load(strcat(rootDir,'Data/Preprocessed/proteomics.mat'))
    % Panc vs All
    per=(LabelIndAll~=2) & (LabelIndAll<5);
    LabelIndAll(per)=1;
    per=(LabelIndAll<5);
    D=LabelIndAll(per);
    D=squareform(pdist(D));
    m=318;
    
    % C=squareform(pdist(A(:,per)'));
    % CorrPermDistTest(C,D,rep*2,'ProteomicsOvarvsAll');
    %ind=181;
    % C=A(:,per)';
    % for i=1:318
    % C(:,i)=C(:,i)/norm(C(:,i));
    % end
    C=squareform(pdist(A(ind,per)'));
    CorrPermDistTest(C,D,rep,'ProteomicsPancvsAllNeuroganin');
    
    % Screening Panc vs All
    pMGC=zeros(m,1);pD=zeros(m,1);pM=zeros(m,1);pP=zeros(m,1);pHHG=zeros(m,1);testMGC=zeros(m,1);testD=zeros(m,1);testM=zeros(m,1);testHHG=zeros(m,1);
    for i=1:m
        C=squareform(pdist(A(i,per)'));
        [pMGC(i),pD(i),pM(i),pP(i), pHHG(i),testMGC(i),testD(i),testM(i),testHHG(i)]=CorrPermDistTest(C,D,rep);
    end
    
    testMGC(:,6)=pMGC;
    testMGC(:,7)=pP;
    testMGC(:,8)=pD;
    testMGC(:,9)=pM;
    testMGC(:,10)=pHHG;
    testMGC(:,2)=testD;
    testMGC(:,3)=testD;
    testMGC(:,4)=testM;
    testMGC(:,5)=testHHG;
    save(strcat(rootDir,'Data/Results/','ScreeningPancvsAll'),'testMGC');
end

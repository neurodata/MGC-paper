function [score,rk]=run_realData_Match(select)
% run real data experiments
% run_realData(10000)  % the replicates for draft
% if nargin<1
%     rep=100;
% end
% if nargin<2
%     select=7;
% end

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
    [score,rk,score2]=MGCGeometry(distC,distP)
%     score2
%     x=cmdscale(distC,1);y=cmdscale(distP,1);   
% [A,B,RX,RY]=MGCDistTransform(distC,distP);
% k=13;
% l=33;
% RX(RX<=k)=1;
% RY(RY<=l)=1;
% RX(RX>k)=0.00001;
% RY(RY>l)=0.00001;
%     x=cmdscale(distC,1);
    %y=mdscale(distP,1, 'Weights',RY);   
%     hold on
%     plot(x,y,'b.');
%     z=log(x.^2);
%     z=z*corr(z,y)/std(z);
%     z=z+mean(y)-mean(z);
%     plot(x,z,'r.');
%     hold off
%     norm(y-z)
%     
% %     z=log(x.^2);
% %     z=z*corr(z,y)/std(z);
% %     z=z+mean(y)-mean(z);
% %     norm(y-z)
% %     z=abs(x).^(0.25);
% %     z=z*corr(z,y)/std(z);
% %     z=z+mean(y)-mean(z);
% %     norm(y-z)
%     z=x;
%     z=z*corr(z,y)/std(z);
%     z=z+mean(y)-mean(z);
%     norm(y-z)
end

%%%
if select==2 || select==7
    load(strcat(rootDir,'Data/Preprocessed/BrainHippoShape.mat'))
    n=114;
    y=squareform(pdist((Label+unifrnd(0,0.01,n,1))));
    [score,rk,score2]=MGCGeometry(LMRS,y)
%     x=cmdscale(LMRS,1);y=cmdscale(y,1);
%     hold on
%     plot(x,y,'b.');
%     z=log(x.^2);
%     plot(x,z,'r.');
%     hold off
    %y=squareform(pdist(Label));
    %y=y+1;
    %for i=1:n
    %    y(i,i)=0;
    %end
    %CorrPermDistTest(LMLS,y,rep,'BrainLMLxY');
%     CorrPermDistTest(LMRS,y,rep, 'BrainLMRxY');
    %CorrPermDistTest(LMLS,LMRS,rep,'BrainLMLxLMR');
end

%%%
if select==3 || select==7
    load(strcat(rootDir,'Data/Preprocessed/semipar.mat'))
    n=109;
    distCCI=squareform(pdist(cci));
    [score,rk,score2]=MGCGeometry(distMigrain(ind,ind),distCCI(ind,ind))
%     x=cmdscale(distMigrain(ind,ind),1);y=cmdscale(distCCI(ind,ind),1);
%     hold on
%     plot(x,y,'b.');
%     z=x;
%     z=z*corr(z,y)/std(z);
%     z=z+mean(y)-mean(z);
%     plot(x,z,'r.');
%     hold off
%     norm(y-z)
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
        %i
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
    ind=181;%296 for hhg
    C=squareform(pdist(A(ind,per)'));
    CorrPermDistTest(C,D,rep,'ProteomicsPancvsNormalNeuroganin');
    
    % Screening Panc vs Normal
    pMGC=zeros(m,1);pD=zeros(m,1);pM=zeros(m,1);pP=zeros(m,1);pHHG=zeros(m,1);testMGC=zeros(m,1);testD=zeros(m,1);testM=zeros(m,1);testHHG=zeros(m,1);
    for i=1:m
        %i
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
    
    C=squareform(pdist(A(:,per)'));
    CorrPermDistTest(C,D,rep,'ProteomicsPancvsAll');
    ind=181;
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
    
    LabelIndAll=LabelIndAll(per)-1;
    n=length(LabelIndAll);
    error=zeros(m,1);k=11;
    for i=1:m
        tmpLabel=zeros(n,1);
        C=squareform(pdist(A(i,per)'));
        [~,~,RX,~]=MGCDistTransform(C,D);
        RX(RX==1)=0;
        RX(RX>k+1)=0;
        RX(RX>0)=1;
        for j=1:n
            tmp=sum(RX(:,j).*LabelIndAll);
            if tmp>floor(k/2)
                tmpLabel(j)=1;
            end
        end
        error(i)=mean(abs(LabelIndAll-tmpLabel));
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

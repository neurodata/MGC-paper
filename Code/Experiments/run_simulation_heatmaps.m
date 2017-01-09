function []=run_simulation_heatmaps(figNumber,repP)
% Used to plot figure 1-8 used in tex. Run like

% total is usually 20.
% pre1 specifies the location to load data.
% pre2 specifies the location to save pictures.

%%% File path searching
if nargin<1
    figNumber='HDHeat';
end
if nargin<2
    repP=100;
end
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
strcat(rootDir,'Code/');
addpath(genpath(strcat(rootDir,'Code/')));
pre1=strcat(rootDir,'Data/Results/'); % The folder to locate data

total=20;
Xmin=zeros(total,1);
Xmax=zeros(total,1);
Ymin=zeros(total,1);
Ymax=zeros(total,1);
K=zeros(total,1);
L=zeros(total,1);


if strcmp(figNumber,'1DHeat')
    nn=60;
    for j=1:total
        filename=strcat(pre1,'CorrIndTestType',num2str(j),'N100Dim1.mat');
        load(filename)
        
        if j<total
            %             phmax=max(max(ph));
            testRepeat=1;
            while testRepeat==1;
                j
                [x,y]=CorrSampleGenerator(j,nn,1,1,1);
                [pMGC,statMGC,pLocalCorr,localCorr,optimalInd]=MGCPermutationTest(squareform(pdist(x)),squareform(pdist(y)),repP,'mcor');
                if pMGC<max(0.05,2/repP) && ((sum(optimalInd==nn^2)==0 && j>5) || (sum(optimalInd==nn^2)==1 && j<=5))
                    testRepeat=0;
                end
            end
            
            %     [~,~,~,optimalInd]=FindLargestRectangles((ph>=phmax-powerThres), [0 0 1]);
            %     optimalInd=find(optimalInd==1);
            [J,I]=ind2sub(size(localCorr),optimalInd);
            Ymin(j)=min(I)-1;
            Ymax(j)=max(I)-1;
            Xmin(j)=min(J)-1;
            Xmax(j)=max(J)-1;
            
            localCorr(Xmin:Xmax,Ymin:Ymax)=localCorr(Xmin:Xmax,Ymin:Ymax)+10;
            statMGC=statMGC+10;
            op2=find(localCorr==statMGC);
            if isempty(op2)
                k=Xmax(j)+1;
                l=Ymax(j)+1;
            else
                [k,l]=ind2sub(size(localCorr),op2);
            end
            K(j)=k;
            L(j)=l;
        end
    end
    save(strcat(pre1,'CorrIndTest1DHeat.mat'),'Xmin','Xmax','Ymin','Ymax','K','L');
else
    figNumber='HDHeat';
    for j=1:total
        filename=strcat(pre1,'CorrIndTestDimType',num2str(j),'N100Dim.mat');
        load(filename)
        thres=0.5;
        ind=find(max(powerMGCM,[],1)>=thres,1,'last');
        if isempty(ind)
            ind=1;
        end
        
        if j<total
            %phmax=max(max(ph));
            testRepeat=1;
            while testRepeat==1;
                j
                [x,y]=CorrSampleGenerator(j,n,dimRange(ind),1,0);
                [pMGC,statMGC,pLocalCorr,localCorr,optimalInd]=MGCPermutationTest(squareform(pdist(x)),squareform(pdist(y)),repP,'mcor');
                if pMGC<max(0.05,2/repP) && ((sum(optimalInd==n^2)==0 && j>5) || (sum(optimalInd==n^2)==1 && j<=5))
                    testRepeat=0;
                end
            end
            
            %     [~,~,~,optimalInd]=FindLargestRectangles((ph>=phmax-powerThres), [0 0 1]);
            %     optimalInd=find(optimalInd==1);
            [J,I]=ind2sub(size(localCorr),optimalInd);
            Ymin(j)=min(I)-1;
            Ymax(j)=max(I)-1;
            Xmin(j)=min(J)-1;
            Xmax(j)=max(J)-1;
            
            localCorr(Xmin:Xmax,Ymin:Ymax)=localCorr(Xmin:Xmax,Ymin:Ymax)+10;
            statMGC=statMGC+10;
            op2=find(localCorr==statMGC);
            if isempty(op2)
                k=Xmax(j)+1;
                l=Ymax(j)+1;
            else
                [k,l]=ind2sub(size(localCorr),op2);
            end
            K(j)=k;
            L(j)=l;
        end
    end
end

save(strcat(pre1,'CorrIndTest',figNumber,'.mat'),'Xmin','Xmax','Ymin','Ymax','K','L');

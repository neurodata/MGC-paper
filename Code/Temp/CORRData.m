function []=CORRData(part,rep,pre1,alpha)

% CORRData(1)

% if nargin<1
%     aList = {'BNU1','BNU2','BNU3','DC1','HNU1','IACAS','IBATRT','IPCAS1','IPCAS2','IPCAS5','IPCAS6','IPCAS8','JHNU','KKI21','LMU3','MPG1','MRN','NKI24mx645','NKI24mx1440','NKI24std2500','NYU1','NYU2','SWU1','SWU2','SWU3','SWU4','UM','UPSM1','Utah1','UWM','XHCUMS'};
% end
if nargin<2
    rep=200; %1000
end
if nargin<3
    pre1='../../../../Data/CORR/';
end
if nargin<4
    alpha=0.05;
end
switch part
    case 1
        aList = {'BNU1','BNU2','BNU3'};
    case 2
        aList =  {'DC1','HNU1','IACAS'};
    case 3
        aList= {'IBATRT','IPCAS1','IPCAS2'};
    case 4
        aList= {'IPCAS5','JHNU','KKI21'};
    case 5
        aList= {'LMU3','MPG1','MRN'};
    case 6
        aList= {'NYU1','NYU2','SWU3','SWU1'};
    case 7
        aList=  {'SWU4','UM','UPSM1','SWU2'};
    case 8
        aList= {'Utah1','UWM','XHCUMS'};
    case 9
        aList= {'IPCAS2','UM','BNU2','HNU1'};
end

option=[0,2,0,0];
total=0;
for l=1:length(aList);
    str1=aList(l);str1=str1{1,1};
    fileDir=strcat(pre1,str1,'/');
    allFiles = dir( fileDir );
    allNames = { allFiles.name };
    indStr = strfind(allNames,'session_1');
    indStr = find(~cellfun(@isempty,indStr));
    allNames=allNames(indStr);
    
    n=length(allNames);
    fileName=allNames(1);
    load(strcat(fileDir,fileName{1,1}));
    [region,timestep]=size(roi_data);
    
    X=zeros(n,timestep,region);
    Y=mvnrnd(0,1,timestep);
    distP=squareform(pdist(Y));
    power=zeros(7,1);
    p1=zeros(7,1);
    total=total+n;
    
    for j=1:n
        fileName=allNames(j);
        load(strcat(fileDir,fileName{1,1}));
        X(j,:,:)=roi_data(:,1:timestep)';
    end
    
    for r=1:region
%         r
        distC=squareform(pdist(X(:,:,r)'));
        %distP=squareform(pdist(X2(:,:,i)));
        [p1(1), p1(2), p1(3), p1(4),p1(5),p1(6),p1(7)]=CorrPermDistTest(distC,distP,rep,'tmp',option);
        for tt=1:7
            if p1(tt)<alpha
                power(tt)=power(tt)+1/region;
            end
        end
    end
    
    fileS=strcat(pre1,str1,'FalseDetection.mat');
    save(fileS,'power','alpha','n','timestep','region','str1');
end
total
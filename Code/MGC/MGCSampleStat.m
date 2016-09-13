function statMGC=MGCSampleStat(localCorr)
% Author: Cencheng Shen
% This function estimate the Oracle MGC (i.e., optimal local correlation)
% from the local correlation map, which we call sample MGC statistic.
%
% It finds a sufficiently large region in the correlation map, such that
% all correlations within the region are significantly large and
% monotonically changing. Then a relatively large correlation in the region
% is used as sample MGC. If no such region exists, use the global
% correlation instead.
%
% Input: a size m*n local correlation map.
% Output: the sample MGC statistic within [-1,1].
[m,n]=size(localCorr);
negCorr=localCorr(2:end,2:end);
negCorr=negCorr(negCorr<0); % negative correlations
eps=3.5*max(norm(negCorr,'fro')/sqrt(length(negCorr)),0.01); % threshold for significantly large correlation

R=(localCorr>eps); % find a region with significantly large correlation
thres=min(2/min(m,n),0.05); % threshold for sufficiently large region
statMGC=localCorr(end); % take global correlation by default

if mean(mean(R))>=thres
    R=Monotone(localCorr,R); % put monotonically changing restriction to reduce the region R
    if mean(mean(R))>0
        [k,l]=find((localCorr>=max(max(localCorr(R))))&(R==1)); % find the scale within R that has the maximum correlation
        ln=ceil(0.1*m); % boundary for checking adjacent rows
        km=ceil(0.1*n); % boundary for checking adjacent columns
        for i=1:length(k)
            ki=k(i);
            li=l(i);
            
            % ensure the adjacent rows does not exceed the local correlation size, same for columns
            left=max(2,li-ln); 
            right=min(n,li+ln);
            upper=max(2,ki-km);
            down=min(m,ki+km);
            
            tmp1=min(localCorr(upper:down,li)); % minimal correlation at given row and adjacent columns
            tmp2=min(localCorr(ki,left:right)); % minimal correlation at given column and adjacent rows
            tmp=max(tmp1,tmp2); % take the max of the minimal correlations, which should be relatively large but does not inflate much over the true optimal statistic
            if tmp>statMGC
                statMGC=tmp; 
            end
        end
    end
end

function [R]=Monotone(localCorr,R)
% An auxiliary function that computes regions in R with monotonically increasing
% or decreasing correlations along the row or column
if nargin<2
    R=ones(size(localCorr));
end
[m,n]=size(localCorr);

R2=cell(4,1);
PD1=zeros(m,n);
PD2=zeros(m,n);
for i=2:m
    PD1(i,2:end)=diff(localCorr(i,:)); % difference within each row
end
for i=2:n
    PD2(2:end,i)=diff(localCorr(:,i)); % difference within each column
end

R2{1}=(PD1>=0)&R; % monotonically increasing along the row within R
R2{2}=(PD1<=0)&R; % monotonically decreasing along the row within R
R2{3}=(PD2>=0)&R; % monotonically increasing along the column within R
R2{4}=(PD2<=0)&R; % monotonically decreasing along the column within R
R=false(m,n);

warning('off','all');
for i=1:4
    t=sum(sum(R2{i}));
    if t>0
        R2{i}=bwareafilt(R2{i},1); % largest connected component of each region 
        R= (R | R2{i}); % combine all monotonically changing and significant regions
    end
end
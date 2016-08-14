function [p,indAll]=MGCScaleVerify(P,thres,alpha)
% An auxiliary function to verify and estimate the MGC optimal scale based
% on the p-values of all local correlations
if nargin<2
thres=0.025;
end
if nargin<3
alpha=0.05;
end

[m,n]=size(P);
%tt=0.025*[m,n];
% tt=[2,2];
% use global p-value by default
p=0.5;
indAll=1;

% If the global p-value is among the top 2.5% of all local p-values, we
% take it directly; otherwise find a rectangular region based on
% monotoneically changing p-values.
if sum(sum(P<P(end)))/(m-1)/(n-1)<thres
    p=P(end);
    indAll=m*n;
else
%     [~, ~, ~, indAll]=FindLargestRectangles((P<=0.05), [0 0 1],[2,2]);
%     if mean(mean(indAll)) >thres
%         p=0.04;
%         indAll=(indAll==1);
%     end
    
    % Regularize the p-values by 0.005
%     PN=floor(P*200)/200; 
    
    R=MonotoneRegion(P,alpha); % Find the largest monotonic rectangular region along the row p-values   
    R2=MonotoneRegion(P',alpha); % Find the largest monotonic rectangular region along the column p-values
    if sum(sum(R2))>sum(sum(R))
        R=R2'; % Take the larger region between row and columns
    end

    % Take the region if and only if its area is larger than 2.5%
    tmp=mean(mean(R));
    if tmp >thres
%         figure
%         imagesc(P(2:end,2:end))
%         figure
%         imagesc(R)
%         p=min(P(R));
        p=prctile(P(R), ceil(thres/tmp*100)); % Take the top 2.5% p-value in the region
        [~, ~, ~, indAll]=FindLargestRectangles((R==1) & (P<=p), [0 0 1]);
        indAll=(indAll==1);
    end
end

function R=MonotoneRegion(P,alpha)
% Find the largest monotonically decreasing or increasing rectangular region along the row p-values
[m,n]=size(P);
% tt=thres*[m,n]; % default minimal size of a rectangular
tt=[2,2];
pdiff=zeros(m,n);
for i=2:m
    tt=P(i,:);
    pdiff(i,2:end)=diff(tt);
%     pdiffDes(i,2:end)=(ttd<=0); % Verify monotonically decreasing
%     pdiffAes(i,2:end)=(ttd>=0); % Verify monotonically increasing
end
% eps=pdiff(2:end,2:end);
% eps=prctile(abs(eps(eps>0)),thres*100);

% for i=2:m
%     tt=P(i,:);
%     ttd=diff(tt);
%     pdiffDes(i,2:end)=(ttd<=0); % Verify monotonically decreasing
%     pdiffAes(i,2:end)=(ttd>=0); % Verify monotonically increasing
% end

R=(pdiff<=0) & (P<alpha);
R2=(pdiff>=0) & (P<alpha);

[~, ~, ~,R] = FindLargestRectangles(R, [0 0 1],tt); % Find the largest rectangular within the monotonically decreasing region
[~, ~, ~, R2] = FindLargestRectangles(R2, [0 0 1],tt); % Find the largest rectangular within the monotonically increasing region
if sum(sum(R2))>sum(sum(R))
    R=R2; % Take the larger region between increasing / decreasing
end

% [~, ~, ~, R] = FindLargestRectangles(R, [0 0 1],tt); % Find the largest rectangular within the monotonically decreasing region
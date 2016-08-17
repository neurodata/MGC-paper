function [p,indAll]=MGCScaleVerify(P,thres,thres2)
% Author: Cencheng Shen
% This function approximates the MGC optimal scale based
% on the p-values of all local correlations
if nargin<2
    thres=0.04; % The threshold is used to: determine if the global p-value is significant enough, and determine if the rectangular area is significant enough
end
if nargin<3
    thres2=0.005; % The threshold to approximate the monotone p-values change
end
[m,n]=size(P);

thres=max(2/min(m,n),thres);
indAll=1;
% default p-value and optimal scale
p=0.5;

% Directly use the global p-value if it is among the top thres% of all local p-values
if sum(sum(P<P(end)))/(m-1)/(n-1)<thres
    p=P(end);
    indAll=m*n;
else
    % otherwise find a rectangular region based on monotoneically changing p-values 
    lim=[2,2]; % The minimal size of a rectangle
    R=SmoothRegion(P,thres2,lim); % Find the largest monotonic rectangular region along the row p-values   
    R2=SmoothRegion(P',thres2,lim); % Find the largest monotonic rectangular region along the column p-values
    if sum(sum(R2))>sum(sum(R))
        R=R2'; % Take the larger region between rows and columns
    end

    % Take the region if and only if its area is larger than the threshold
    tmp=mean(mean(R(2:end,2:end)));
    if tmp >thres
        p=prctile(P(R), ceil(thres/tmp*100)); % Take the top p-value in the region
        [~, ~, ~, indAll]=FindLargestRectangles((R==1) & (P<=p), [0 0 1],lim);
        indAll=find(indAll==1);
    end
end

function R=SmoothRegion(P,thres,lim)
% Find the largest monotonically decreasing or increasing rectangular region along the row p-values
[m,n]=size(P);
pdiff=zeros(m,n);
for i=2:m
    tt=P(i,:);
    pdiff(i,2:end)=diff(tt);
end

R=(pdiff<thres); % Check monotone decreasing, but also allows small p-value increase no more than the threshold
R2=(pdiff>-thres); % Check monotone increasing, but also allows small p-value decrease no more than the threshold

[~, ~, ~,R] = FindLargestRectangles(R, [0 0 1],lim); % Find the largest rectangular within the (approximately) monotonically decreasing region
[~, ~, ~, R2] = FindLargestRectangles(R2, [0 0 1],lim); % Find the largest rectangular within the (approximately) monotonically increasing region
if sum(sum(R2))>sum(sum(R))
    R=R2; % Take the larger region
end
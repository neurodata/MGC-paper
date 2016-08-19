function [p,indAll]=MGCScaleVerify(P,thres,thres2)
% Author: Cencheng Shen
% This function approximates the p-value and the optimal scale for sample MGC, 
% based on the p-values of all local correlations. 
% The algorithm first looks for a smooth region in the p-value map.
% Then it takes the global p-value directly f it is small enough among all local p-values, 
% otherwise it checks whether the smooth region is large enough and approximate the MGC p-value from the smooth region instead.
% At last, the smooth rectangle that is bounded by the sample MGC p-value is taken as the optimal scales. 
if nargin<2
    thres=0.05; % the threshold is used to: determine if the global p-value is significant enough, and determine if the rectangular area is significant enough
end
if nargin<3
    thres2=0.005; % the threshold to approximate the monotone p-values change
end
[m,n]=size(P);
sz=[2,2]; % minimal rectangle size
thres=max(1/min(m,n),thres); % in case sample size is too small, increase threshold

% default p-value
p=1;
indAll=1;

R=SmoothRegion(P,thres2); % find the largest smooth region in the p-value map
% directly use the global p-value if it is among the top thres% of all local p-values
if sum(sum(P<P(end)))/(m-1)/(n-1)<thres
    p=P(end);
    indAll=m*n;
else
    tmp=mean(mean(R(2:end,2:end)));
    % approximate MGC p-value from the smooth region if and only if the region area is larger than the threshold
    if tmp>thres
        % select a top p-value in the region based on thres and the area of the region.
        p=prctile(P(R), min(ceil(thres/tmp*100),100));
    end
end

% the largest rectangle within the smooth region bounded by the
% p-value is taken as the optimal scale
[~, ~, ~, R]=FindLargestRectangles((R==1) & (P<=p), [0 0 1],sz);
if sum(sum(R(2:end,2:end)))>0
    indAll=find(R==1);
end

function R=SmoothRegion(P,thres)
% Find the largest smooth rectangular region in the p-value map, by considering the
% largest monotonically decreasing or increasing rectangular region along
% the row or column p-values, but allowing small p-value increase or
% decrease as specified by thres. Then further require all p-values in the smooth region
% to be less than 0.5.
[m,n]=size(P);
sz=[2,2]; % minimal rectangle size

PD=cell(2,1);
PD{1}=zeros(m,n); % store the difference within rows
PD{2}=zeros(m,n); % store the difference within columns
for i=2:m
    tt=P(i,:);
    PD{1}(i,2:end)=diff(tt);
end
for i=2:n
    tt=P(:,i);
    PD{2}(2:end,i)=diff(tt);
end

RC=cell(4,1);
R=zeros(m,n);

RC{1}=(PD{1}<thres); % check monotone decreasing , but also allows small p-value increase no more than the threshold
RC{2}=(PD{1}>-thres); % check monotone increasing , but also allows small p-value decrease no more than the threshold
RC{3}=(PD{2}<thres); % repeat for column changes
RC{4}=(PD{2}>-thres); 
                
for i=1:4
    [~, ~, ~,RC{i}] = FindLargestRectangles(RC{i}, [0 0 1],sz); % find the largest rectangles within the (approximately) monotonically decreasing / increasing region
    R=(R==1) | (RC{i}==1); % combine the four rectangular regions
end

[~,~,~,R]=FindLargestRectangles((R==1) & (P<0.5), [0 0 1],sz); % in the smooth region, find the largest rectangle with all p-values less than 0.5
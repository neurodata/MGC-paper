function [p,indAll]=MGCScaleVerify(P,gamma,tau)
% Author: Cencheng Shen
% This function approximates the p-value and the optimal scale for sample MGC, 
% based on the p-value map of all local correlations. 
% The algorithm first looks for a smooth region in the p-value map, and
% takes the largest p-value in the smooth region.
% However, if the global p-value is small enough among all local p-values, 
% the global scales is used for MGC;
% and if there exists a smooth rectangular region that is large enough, 
% a smaller p-value within the smooth rectangle is approximated for MGC.
% Once we determine the sample MGC p-value, the smooth rectangle that is bounded 
% by the p-value is taken as the optimal scales (and further include the global scale if necessary).
% If the optimal scale is empty, we take the largest p-value in the p-value
% map and all scales (from k,l=2 til n) for MGC.
if nargin<2
    gamma=0.1; % gamma is used to: determine if the global p-value is significant enough, determine if the rectangular region is significant enough, and approximate a small p-value in the significant rectangular region
end
if nargin<3
    tau=0.005; % tau is a threshold to approximate the monotone p-values change
end
[m,n]=size(P);
gamma=max(2/min(m,n),gamma); % increase gamma accordingly in case the sample size is too small
indAll=[]; % optimal scale

% find the largest smooth region in the p-value map
R=SmoothRegion(P,tau); 
% p=max(max(P(2:end,2:end)));
if sum(sum(R))>0
    p=max(P(R)); % take the largest p-value in the smooth region first, as long as the smooth region is not empty
end

% further check for global p-value and smooth rectangle
if sum(sum(P<P(end)))/(m-1)/(n-1)<gamma
    p=P(end); % directly use the global p-value if it is among the top 10% of all local p-values
    indAll=m*n;
else
    [~,~,~,R]=FindLargestRectangles(R, [0 0 1],[2,2]); % find the largest rectangle within the smooth region
    tmp=mean(mean(R));
    % approximate a smaller p-value from the smooth rectangle if and only if the area is larger than gamma
    if tmp>gamma
        p=prctile(P(R), gamma/tmp*50); % take a small p-value that is 5/area(R)% of all p-values in the smooth rectangle
    end
end

% the largest smooth rectangle bounded by the sample MGC p-value is taken as the optimal scale
[~, ~, ~, R]=FindLargestRectangles((R==1) & (P<=p), [0 0 1],[2,2]);
indAll=[find(R==1);indAll]; % include the global scale if necessary

% use the largest of all p-values when the optimal scale is empty
if isempty(indAll)
    p=max(max(P(2:end,2:end)));
    indAll=find(P<=p);
end

function R=SmoothRegion(P,tau)
% Find the largest smooth rectangular region in the p-value map, by considering the
% largest monotonically decreasing or increasing scales along
% the row or column p-values, but allowing small p-value increase or
% decrease as specified by tau. 
[m,n]=size(P);

PD=cell(2,1);
PD{1}=zeros(m,n); % store the p-value changes within rows
PD{2}=zeros(m,n); % store the p-value changes within columns
for i=2:m
    tt=P(i,:);
    PD{1}(i,2:end)=diff(tt);
end
for i=2:n
    tt=P(:,i);
    PD{2}(2:end,i)=diff(tt);
end

RC=cell(4,1);

RC{1}=(PD{1}<tau); % check monotone decreasing , but also allows small p-value increase no more than tau
RC{2}=(PD{1}>-tau); % check monotone increasing , but also allows small p-value decrease no more than tau
RC{3}=(PD{2}<tau); % repeat for column changes
RC{4}=(PD{2}>-tau); 
                
for i=1:4
    [~, ~, ~,RC{i}] = FindLargestRectangles(RC{i} & (P<1), [0 0 1],[2,2]); % find the largest rectangle within each (approximately) monotonically decreasing / increasing region
end

R=(RC{1} | RC{2} | RC{3} | RC{4}); % combine all four rectangles together into the smooth region

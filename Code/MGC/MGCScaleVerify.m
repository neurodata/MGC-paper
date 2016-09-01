function [p,indAll]=MGCScaleVerify(P,rep)
% Author: Cencheng Shen
% This function approximates the p-value and the optimal scale for sample MGC,
% based on the p-value map of all local correlations.
% The algorithm first searches for large smooth regions in the p-value map.
% If the global p-value is small enough among all local p-values, the global p-value is used for MGC;
% or if there exists a smooth rectangular region that is large enough,
% a smaller p-value within the smooth rectangle is approximated for MGC.
% Once we determine the sample MGC p-value, the smooth rectangle that is
% bounded above by the p-value is taken as the optimal scales, which shall include the global scale if necessary.

gamma=0.04; % gamma is used to: determine if the global p-value is significant enough, and determine if the rectangular region is significant enough
[m,n]=size(P);
gamma=max(3/min(m,n),gamma); % increase gamma accordingly in case the sample size is too small

% find large smooth regions in the p-value map
R=SmoothRegions(P,rep);

% check for global p-value and smooth rectangle
if sum(sum(P<P(end)))/(m-1)/(n-1)<=gamma
    p=P(end); % directly use the global p-value if it is among the top 100*gamma% of all local p-values
else
    [~,~,~,R]=FindLargestRectangles(R, [0 0 1],[2,2]); % find the largest rectangle within the smooth regions
    tmp=mean(mean(R));
    if tmp>=gamma % if the smooth region is too small, take all scales
        [h,edges]=histcounts(P(R),min(ceil(sum(sum(R))/20),100)); % partition the smooth p-values into bins
        [~,ii]=max(h); % find the bin with most counts
        p=edges(ii+1); % use the largest p-value in that bin for MGC
    else
         p=P(end);
    end
end

% the largest smooth rectangle bounded by the sample MGC p-value is taken as the optimal scale
[~, ~, ~, R]=FindLargestRectangles((R==1) & (P<=p), [0 0 1],[2,2]);
if sum(sum(R))==0
    p=P(end); % if empty, use the global p-value
end

% include the global scale if necessary
indAll=find(R==1);
if (p>=P(end) && R(m*n)==0) 
    indAll=[indAll;m*n];
end

function R=SmoothRegions(P,rep)
% Find the smooth regions in the p-value map, by considering the
% largest monotonically decreasing or increasing scales along
% the row or column p-values, but allowing small p-value increase or
% decrease as specified by tau. 
[m,n]=size(P);
tau=min(2/rep,0.005);

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

RC{1}=(PD{1}<=tau); % check monotone decreasing , but also allows small p-value increase no more than tau
RC{2}=(PD{1}>=-tau); % check monotone increasing , but also allows small p-value decrease no more than tau
RC{3}=(PD{2}<=tau); % repeat for column changes
RC{4}=(PD{2}>=-tau); 
                
for i=1:4
    [~, ~, ~,RC{i}] = FindLargestRectangles(RC{i} & (P<1), [0 0 1],[2,2]); % find the largest rectangle within each (approximately) monotonically decreasing / increasing region
end

R=(RC{1} | RC{2} | RC{3} | RC{4}); % combine all four rectangles together into the smooth regions
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
% Otherwise, we use the mean p-value of all local p-values instead.

gamma=0.05; % gamma is used to: determine if the global p-value is significant enough, determine if the rectangular region is significant enough, and approximate a small p-value in the significant rectangular region
tau=10;
% if nargin<3
%     tau=0.0002; % tau is a threshold to approximate the monotone p-values change for smooth regions
% end
[m,n]=size(P);
%gamma=max(3/min(m,n),gamma); % increase gamma accordingly in case the sample size is too small

% find large smooth regions in the p-value map
R=SmoothRegions(P,rep);
% [h,edges]=histcounts(P(R),ceil(sum(sum(R))/10));
% [~,ii]=max(h);
% p=edges(ii+1);

% check for global p-value and smooth rectangle
if sum(sum(P<P(end)))/(m-1)/(n-1)<=gamma
    p=P(end); % directly use the global p-value if it is among the top 100*gamma% of all local p-values
else
     [~,~,~,R]=FindLargestRectangles(R, [0 0 1],[2,2]); % find the largest rectangle within the smooth regions
%     tmp=mean(mean(R));
    % approximate a smaller p-value from the smooth rectangle if and only if the area is no smaller than gamma
    %if tmp>=max(2/min(m,n),gamma) % increase gamma accordingly in case the sample size is too small
        % take a small p-value that is 100*gamma/2area(R)% of all p-values in the smooth rectangle.
        % For example, if the area of R equals gamma, the median
        % p-value within R is use for the MGC p-value.
        %         p=prctile(P(R), gamma/tmp*100);
        [h,edges]=histcounts(P(R),min(ceil(sum(sum(R))/tau),100));
        [~,ii]=max(h);
%         if h(ii)>20
           p=edges(ii+1);
%     else
%         p=mean(P(P<1)); % otherwise, use the mean p-value for MGC
%     end
end

% the largest smooth rectangle bounded by the sample MGC p-value is taken as the optimal scale
[~, ~, ~, R]=FindLargestRectangles((R==1) & (P<=p), [0 0 1]);
% if mean(mean(R))<0.01 && p~=P(end)
if sum(sum(R))<max(m*n*gamma/5,tau) && p~=P(end)
    % if the scales from above return empty, relax the smooth region contraint for optimal scale. 
    % This rarely happens except when the mean p-value is taken.
%     p=P(end);
%    p=mean(P(P<1)); 
    [~, ~, ~, R]=FindLargestRectangles((P<=p), [0 0 1]);
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
% tau=[PD{1}(PD{1}>-1);PD{1}(PD{1}>-1)];
% tau=prctile(abs(tau),10)
% tau=min(tau,0.005);

RC{1}=(PD{1}<=tau); % check monotone decreasing , but also allows small p-value increase no more than tau
RC{2}=(PD{1}>=-tau); % check monotone increasing , but also allows small p-value decrease no more than tau
RC{3}=(PD{2}<=tau); % repeat for column changes
RC{4}=(PD{2}>=-tau); 
                
for i=1:4
    [~, ~, ~,RC{i}] = FindLargestRectangles(RC{i} & (P<1), [0 0 1]); % find the largest rectangle within each (approximately) monotonically decreasing / increasing region
end

R=(RC{1} | RC{2} | RC{3} | RC{4}); % combine all four rectangles together into the smooth regions
function [p,indAll]=MGCScaleVerify(P,thres,thres2)
% Author: Cencheng Shen
% This function approximates the optimal scale and p-value of MGC based
% on the p-values of all local correlations. 
% The global p-value is used directly if it is small enough among all local
% p-values, otherwise the algorithm looks for a sufficiently large and smooth region in the
% p-value map and approximate the MGC p-value from the smooth region.
% The smooth rectangle that is bounded by the MGC p-value is taken as the optimal scales. 
if nargin<2
    thres=0.05; % The threshold is used to: determine if the global p-value is significant enough, and determine if the rectangular area is significant enough
end
if nargin<3
    thres2=0.005; % The threshold to approximate the monotone p-values change
end
[m,n]=size(P);
lim=[2,2]; % The minimal size of a rectangle

thres=max(1/min(m,n),thres);

% default p-value and optimal scale
p=0.5;
indAll=1;

warning('off','all');
% Directly use the global p-value if it is among the top thres% of all local p-values
if sum(sum(P<P(end)))/(m-1)/(n-1)<thres
    p=P(end);
    % Identify optimal scales
    [~, ~, ~, R]=FindLargestRectangles((P<=p), [0 0 1],lim);
    if sum(sum(R(2:end,2:end)))>=lim(1)*lim(2)
        indAll=find(R==1);
    else
        indAll=m*n;
    end
else
    R=SmoothRegion(P,thres2,lim); % Find the largest smooth region in the p-value map
    
    tmp=mean(mean(R(2:end,2:end)));
    % Take the smooth region if and only if the area is larger than the threshold
    if tmp>thres
        % Select a top p-value in the region based on thres and the area of the region.
        % Alternatively, one can go through a list of candidate p-values, or the significance level alpha
        p=prctile(P(R), min(ceil(thres/tmp*100),100));

        % Find the largest rectangle within the smooth region that are
        % covered by the candidate p-value
        [~, ~, ~, R]=FindLargestRectangles((R==1) & (P<=p), [0 0 1],lim);
        indAll=find(R==1);
%         if p<0.05
%             tmp
%             p
%             P(end)
%         end
    end
end
%                     figure
%                     imagesc(R)

function R=SmoothRegion(P,thres,lim)
% Find the largest smooth rectangular region in the p-value map, by considering the
% largest monotonically decreasing or increasing rectangular region along
% the row or column p-values, but allowing small p-value increase or
% decrease as specified by thres. Then further require all p-values in the smooth region
% to be less than 0.5.
[m,n]=size(P);
PD=cell(2,1);
PD{1}=zeros(m,n); % Store the difference within rows
PD{2}=zeros(m,n); % Store the difference within columns
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

RC{1}=(PD{1}<thres); % Check monotone decreasing , but also allows small p-value increase no more than the threshold
RC{2}=(PD{1}>-thres); % Check monotone increasing , but also allows small p-value decrease no more than the threshold
RC{3}=(PD{2}<thres); % Repeat for column changes
RC{4}=(PD{2}>-thres); 
                
for i=1:4
    [~, ~, ~,RC{i}] = FindLargestRectangles(RC{i}, [0 0 1],lim); % Find the largest rectangles within the (approximately) monotonically decreasing / increasing region
    R=(R==1) | (RC{i}==1); % Combine the four rectangular regions
end

[~,~,~,R]=FindLargestRectangles((R==1) & (P<0.5), [0 0 1],lim); % In the smooth region, find the largest rectangle with all p-values less than 0.5
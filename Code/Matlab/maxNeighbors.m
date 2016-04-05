function [n1]=maxNeighbors(power1,all,dCor1N,dCor1A)
% return the index corresponding to the maximal value in power1;
% if there are more than one such indices, return the first such scale that maximizes
% the mean difference between the the test statistics under the alternative
% and the null.
if nargin<2
    all=1; % Set all to 1 to return all indices in case of multiple maximal values; otherwise the last index is returned
end
if nargin<3
    dCor1N=0;
    dCor1A=0;
    testMean=false;
else
    testMean=true; 
end
pmax=max(max(power1)); % Maximal power
ind=find(power1==pmax);
n1=ind;

% Check if there are multiple indices that maximize the testing power
if all==0
    if testMean==true
        meanInd=mean(dCor1A,3)-mean(dCor1N,3);
        meanInd=meanInd(ind);
        pmax=max(meanInd);
        ind2=find(meanInd==pmax);
        if (isempty(ind2)==false)
            ind2=ind2(1);
            ind=ind(ind2);
        end
    end
    n1=ind(end);
end
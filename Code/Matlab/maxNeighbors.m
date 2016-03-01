function [n1,k,l]=maxNeighbors(power1,dCor1N,dCor1A)
% return the index corresponding to the maximal value in power1;
% if there are more than one such indices, return the first such scale that maximizes
% the mean difference between the the test statistics under the alternative
% and the null.
if nargin<2
    dCor1N=0;
    dCor1A=0;
    testMean=false;
else
    testMean=true;
end
n=size(power1,1);
pmax=max(max(power1)); % Maximal power
ind=find(power1==pmax);
% n1=ind(1);

% if length(ind)>1
%     l=ceil(n1/n);k=n1-n*(l-1);
%     pmax=k*l;
%     for i=2:length(ind);
%         nt=ind(i);
%         l=ceil(nt/n);k=nt-n*(l-1);
%         if k*l>=pmax;
%             pmax=k*l;
%             n1=ind(i);
%         end
%     end
% end

% % Check if there are multiple indices that maximize the testing power
if length(ind)>1 && testMean==true
    meanInd=mean(dCor1A,3)-mean(dCor1N,3);
    meanInd=meanInd(ind);
    pmax=max(meanInd);
    ind2=find(meanInd==pmax);
    if (isempty(ind2)==false)
        ind2=ind2(1);
        ind=ind(ind2);
    end
end

% n1 is the linear indexing of the matrix power1,
% [k,l] is the row and column subscript of n1.
n1=ind(1);
l=ceil(n1/n);k=n1-n*(l-1);
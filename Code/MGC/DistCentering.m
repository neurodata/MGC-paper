function [A]=DistCentering(X,option)
% An auxiliary function that properly centers the distance matrix X,
% depending on the choice of global corr.
n=size(X,1);

% centering for dcorr / mcorr / Mantel
switch option
    case 'dcor'
        EX=repmat(mean(X,1),n,1); % column centering
        %%% EX=repmat(mean(X,1),n,1)+repmat(mean(X,2),1,n)-mean(mean(X)); % this is original double-centering
    case 'mcor'
        EX=repmat(sum(X,1)/n,n,1);
        EX=EX+X/n;
    case 'mantel'
        EX=sum(sum(X))/n/(n-1);
end
A=X-EX;

% for mcorr or Mantel, exclude the diagonal entries
if strcmp(option,'mcor') || strcmp(option,'mantel')
    %%% meanX=sum(sum(X))/n^2;
    for j=1:n
        A(j,j)=0;
        %%% A(j,j)=sqrt(2/(n-2))*(mean(X(:,j))-meanX)*1i;  % the original diagonal modification of mcorr
    end
end
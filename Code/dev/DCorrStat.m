function dcorr=DCorrStat(C,D)

C=centering(C);
D=centering(D);
dcorr=sum(sum(C.*D))/sqrt(sum(sum(C.*C))*sum(sum(D.*D)));

function A=centering(X)
[n,m]=size(X);
EX=repmat(mean(X,1),n,1)+repmat(mean(X,2),1,n)-mean(mean(X));
EX=EX+X/n;
A=X-EX;
for j=1:m
    A(j,j)=0;
    %%% A(j,j)=sqrt(2/(n-2))*(mean(X(:,j))-meanX)*1i;  % the original diagonal modification of mcorr
end
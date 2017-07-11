test=zeros(21,1);
for i=1:21
d=1;
n=100;
rho=(i-1)/20;
cov1=[eye(d) rho*ones(d)];
cov2=[rho*ones(d) eye(d)];
covT=[cov1' cov2'];
x=mvnrnd(zeros(n,2*d),covT,n);
y=x(:,d+1:2*d);
x=x(:,1:d);
A=squareform(pdist(x));
B=squareform(pdist(y));
tmp=MGCLocalCorr(A,B,'mcor');
test(i)=tmp(end);
end
function [testStat] = HSIC(x,y)
%If this is 1, save the test threshold for later use
% doSaveThresh = 0;
K=HSIC_Kernel(x);
L=HSIC_Kernel(y);
m=size(K,1);
% dx=size(X,2);

H = eye(m)-1/m*ones(m,m);

Kc = H*K*H;
testStat = 1/m^2*sum(sum(Kc'.*L));

function [H]=HSIC_Kernel(X)

% H=X*X';
size1=size(X,1);
if size1>100
    Xmed = X(1:100,:);
    size1 = 100;
else
    Xmed = X;
end
G = sum((Xmed.*Xmed),2);
Q = repmat(G,1,size1);
R = repmat(G',size1,1);
dists = Q + R - 2*Xmed*Xmed';
dists = dists-tril(dists);
dists=reshape(dists,size1^2,1);
deg = sqrt(0.5*median(dists(dists>0)));  %rbf_dot has factor of two in kernel

%Note : patterns are transposed for compatibility with C code.
patterns1=X;
patterns2=X;
size1=size(patterns1);
size2=size(patterns2);

G = sum((patterns1.*patterns1),2);
H = sum((patterns2.*patterns2),2);

Q = repmat(G,1,size2(1));
R = repmat(H',size1(1),1);

H = Q + R - 2*patterns1*patterns2';

H=exp(-H/2/deg^2);

% HI=sum(H,2);
% H=H./repmat(HI,1,size(H,2));
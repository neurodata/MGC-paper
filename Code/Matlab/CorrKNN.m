function error=CorrKNN(C,Label,k)

rep1=1000;
n=size(C,1);
y=squareform(pdist(Label));
%y=(y>0)+1;
y=y+1;
for i=1:n
    y(i,i)=0;
end
error=zeros(k,1);

% [p1,p2,p3,~,~,~,~,neighbor]=CorrPermDistTest(C,y,rep1,0);
% l=ceil(neighbor(1)/n);
% k=neighbor(1)-(l-1)*n;
for i=1:n
    ind=[1:i-1,i+1:n];
    cc=C(ind,i);
    [~,indT]=sort(cc,'ascend');
    for j=1:k
        ind=indT(1:j);
        testLabel=mode(Label(ind));
        if testLabel~=Label(i);
            error(j)=error(j)+1/n;
        end
    end
end

plot(1:k,error,'.-');
    
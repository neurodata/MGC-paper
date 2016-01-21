clear
power=zeros(7,1);
p=zeros(7,1);
rep1=100;rep2=100;

str1='..\..\Data\BNU1\BNU1_00';
str2=25863;
str3='_session_1_timeseries.npz.mat';
ts=200;
n=50;

% str1='..\..\Data\BNU2\BNU2_00';
% str2=25920;
% str3='_session_1_timeseries.npz.mat';
% ts=240;
% n=50;

alpha=0.05;
roin=200;
X=zeros(n,ts,roin);
Y=mvnrnd(0,1,n);
distP=squareform(pdist(Y));

for i=1:n
    filename=strcat(str1,num2str(str2+i),str3);
    load(filename);
    X(i,:,:)=roi_data';
end

for i=1:roin
    distC=squareform(pdist(X(:,:,i)));
    type=[distC distP];
    [p(1), p(2), p(3), p(4),p(5),p(6),p(7)]=CorrPermDistTest(distC,distP,rep1,rep2);
    for j=1:7
        if p(j)<alpha
            power(j)=power(j)+1;
        end
    end
end
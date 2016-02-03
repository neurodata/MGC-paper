clear
power=zeros(7,1);
p1=zeros(7,1);
rep1=1000;rep2=1000;

str1='..\..\Data\BNU1\BNU1_00';
str2=25863;
str3='_session_1_timeseries.npz.mat';
str4='_session_2_timeseries.npz.mat';
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
X2=zeros(n,ts,roin);
Y=mvnrnd(0,1,ts);
distP=squareform(pdist(Y));

for i=1:n
    filename=strcat(str1,num2str(str2+i),str3);
    load(filename);
    X(i,:,:)=roi_data';
    filename=strcat(str1,num2str(str2+i),str4);
    load(filename);
    X2(i,:,:)=roi_data';
end

% %%% Estimate optimal scales jointly All
% power1=0;power2=0;power3=0;
% for i=1:roin
%     distC1=squareform(pdist(X(:,:,i)));
%     [p1,p2,p3]=CorrPermDistTest(distC1,distP,rep1,0,'BNU1');
%     power1=power1+p1/roin;power2=power2+p2/roin;power3=power3+p3/roin;
% end
% neighbor=zeros(3,1);
% neighbor(1)=verifyNeighbors(1-power1);neighbor(2)=verifyNeighbors(1-power2);neighbor(3)=verifyNeighbors(1-power3);
% 
% p1=zeros(7,1);
% for i=1:roin
%     distC1=squareform(pdist(X(:,:,i)));
%     [p1(1), p1(2), p1(3), p1(4),p1(5),p1(6),p1(7)]=CorrPermDistTest(distC1,distP,0,rep2,'BNU11',neighbor);
%     for j=1:7
%         if p1(j)<alpha
%             power(j)=power(j)+1/roin;
%         end
%     end
% end

%%% Estimate optimal scales jointly 2
for i=1:roin
    distC1=squareform(pdist(X(:,:,i)'));
    distC2=squareform(pdist(X2(:,:,i)'));
    [power1,power2,power3]=CorrPermDistTest(distC1,distP,rep1,0,'BNU1');
    %[power4,power5,power6]=CorrPermDistTest(distC2,distP,rep1,0,'BNU1');
    neighbor=zeros(3,1);
    neighbor(1)=verifyNeighbors(1-power1);neighbor(2)=verifyNeighbors(1-power2);neighbor(3)=verifyNeighbors(1-power3);
%     [p1(1), p1(2), p1(3), p1(4),p1(5),p1(6),p1(7)]=CorrPermDistTest(distC1,distP,0,rep2,'BNU11',neighbor);
%     for j=1:7
%         if p1(j)<alpha
%             power(j)=power(j)+1/2/roin;
%         end
%     end
    [p1(1), p1(2), p1(3), p1(4),p1(5),p1(6),p1(7)]=CorrPermDistTest(distC2,distP,0,rep2,'BNU12',neighbor);
    for j=1:7
        if p1(j)<alpha
            power(j)=power(j)+1/roin;
        end
    end
end

% %%% Estimate optimal scales separately
% for i=1:roin
%     distC1=squareform(pdist(X(:,:,i)));
%     [p1(1), p1(2), p1(3), p1(4),p1(5),p1(6),p1(7)]=CorrPermDistTest(distC1,distP,rep1,rep2,'BNU10');
%     for j=1:7
%         if p1(j)<alpha
%             power(j)=power(j)+1/2/roin;
%         end
%     end
%     
%     distC2=squareform(pdist(X2(:,:,i)));
%     [p1(1), p1(2), p1(3), p1(4),p1(5),p1(6),p1(7)]=CorrPermDistTest(distC2,distP,rep1,rep2,'BNU10');
%     for j=1:7
%         if p1(j)<alpha
%             power(j)=power(j)+1/2/roin;
%         end
%     end
% end
% 
% 
%%% Estimate optimal scales separately
for i=1:roin
    distC=squareform(pdist(X(:,:,i)'));
    %distP=squareform(pdist(X2(:,:,i)));
    [p1(1), p1(2), p1(3), p1(4),p1(5),p1(6),p1(7)]=CorrPermDistTest(distC,distP,rep1,rep2);
    for j=1:7
        if p1(j)<alpha
            power(j)=power(j)+1/roin;
        end
    end
end
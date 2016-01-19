clear
ts=200;
alpha=0.05;
rep1=200;rep2=200;
distP=squareform(pdist(mvnrnd(0,1,ts)));
power=zeros(7,1);
p=zeros(7,1);

str1='..\..\Data\BNU1\BNU1_00';
str2=25863;
str3='_session_1_timeseries.npz.mat';
lim=51;

for i=1:lim
    filename=strcat(str1,num2str(str2+i),str3);
    load(filename);
    distC=squareform(pdist(roi_data'));
    type=[distC distP];
    [p(1), p(2), p(3), p(4),p(5),p(6),p(7)]=CorrPermDistTest([distC distP],rep1,rep2);
    for j=1:7
        if p(j)<alpha
            power(j)=power(j)+1;
        end
    end
end
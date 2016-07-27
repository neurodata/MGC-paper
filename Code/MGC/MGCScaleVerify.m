function [p,indAll,t]=MGCScaleVerify(P,alpha)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations
if nargin<2
    alpha=0.05;
end
[p1]=Validate(P,alpha);
[p2]=Validate(P',alpha);
p=min(p1,p2);
indAll=find(P<=p);
if (p<=0.05 && P(end)>0.05)
    p
    P(end)
%     figure
%     imagesc(P);
%     caxis([0 0.05])
end

function [p]=Validate(P,alpha)
alpha=0.05;
[m,n]=size(P);
nn=ceil(0.05*n);
%  nn=1;
% Find motonocally decreasing p-values on each column
pdiff=zeros(m-2,n);
for i=2:n
    tt=P(2:end,i);
    ttd=diff(tt);
    tt=tt(2:end);
    if  (min(tt)<alpha)
          pdiff(:,i)=(ttd<=0); %& (tt<alpha);
          %pdiff(:,i-1)=(tt<alpha);
    end
end
% pdiff=[ones(m-2,nn) pdiff];
% pdiff=[pdiff ones(m-2,nn)];
%
%nn=1;
pp=zeros(n,1);
sind=zeros(n,2);
for i=2:n
    is=max(2,i-nn);
    ie=min(n,i+nn);
    tmp=sum(pdiff(:,is:ie),2);
    tmp=find(tmp<ie-is+1);
    tmp=[1;tmp;m+1];
    tmp2=diff(tmp);
    %ind=find(tmp2==max(tmp2),1,'last');
    jm=max(tmp2);
    j=find(tmp2==jm,1,'last');
    
%     if sum(P(tmp(j+1)-1,i-1:i+1)<alpha)==3
        pp(i)=(jm-1)/(m-2);
        sind(i,:)=[tmp(j)+1, tmp(j+1)-1];
%     end
%      else
%           sind(i,:)=[1,1];
%      end
end

[mm,i]=max(pp);
% i
% figure
% plot(pp);
%thres=0.4;
thres=max(0.2, pp(n)+2*std(pp(2:n)));
% mm
% thres
if mm>=thres;
%     mm
%     thres
    i=find(pp>=thres);
    p=1;
    for t=1:length(i)
        j=i(t);
         pt=min(min(P(sind(j,1):sind(j,2),j)));
         p=min(p,pt);
    end
else
    p=min(P(end),median(P(P<1)));
end
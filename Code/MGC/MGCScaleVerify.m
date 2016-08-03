function [p,indAll,t]=MGCScaleVerify(P,alpha)
% % An auxiliary function to verify and estimate the MGC optimal scale based
% % on the p-values of all local correlations
if nargin<2
    alpha=0.05;
end

% end
[m,n]=size(P);
p=P(end);

indAll=m*n;
pc=mean(mean(P<p))*100;
if pc<=5
    prc=pc;
else
    prc=5:5:50;
    prc=[prc';pc;mean(mean(P<(alpha-0.001)))*100];
    prc=sort(prc,'ascend');
end
for j=1:length(prc)
    alpha=prctile(P(P<1), prc(j));
    try
        CC=bwareafilt(P<alpha,1);
        CC=ValidateRow(CC);
%         CC=ValidateRow(CC');
        
        CCC=ValidateRow(CC');
        if sum(sum(CCC))>sum(sum(CC))
            CC=CCC;
        end
        CC=bwareafilt(CC,1);    
    catch
        CC=zeros(m,n);
    end
    mm=sum(sum(CC))/(m-1)/(n-1);
    thres=max([prc(j)*0.3/100;0.05]);
    if mm>thres
        indAll=(CC==1);
        p=alpha;
        break;
    end
end

% [p1]=Validate(P,rep);
% [p2]=Validate(P',rep);
% p=min(p1,p2);
% indAll=find(P<=p);
if (p<0.05 && P(end)>0.05)
%     mm
%     thres
    %     sum(sum(CC))/(m-1)/(n-1)
%     figure
%     imagesc(CC1)
%     figure
%     imagesc(CC')
    p
    P(end)
end

function CC2=ValidateRow(CC)
[m,n]=size(CC);
nn=ceil(0.02*m);
% nn=1;
CC2=CC;
for i=2:m
    is=max(2,i-nn);
    ie=min(m,i+nn);
    sind=zeros(1,n);
    tmp=find(sum(CC(is:ie,:),1)<ie-is+1);
    tmp=[1;tmp';n+1];
    tmp2=diff(tmp);
    j=find(tmp2==max(tmp2),1,'last');
    %if max(tmp2)>ceil(0.02*n)
    if max(tmp2)>ceil(0.1*n)
       sind(tmp(j)+1:tmp(j+1)-1)=1;
    end
    CC2(i,:)=sind;
end

function [p]=Validate(P,rep)
prc=5:5:50;
[m,n]=size(P);
p=P(end);
for j=1:length(prc)
    alpha=prctile(P(P<1), prc(j));
%     if alpha*j>p
%         break;
%     end
   nn=ceil(0.025*n);
%     nn=ceil(0.05*n);
%        nn=1;
     
     pdiff=zeros(m-2,n);
     for i=2:n
         tt=P(2:end,i);
           ttd=diff(tt);
              tt=tt(2:end);
%           if  (min(tt)<alpha)
              pdiff(:,i)=(ttd<=0) & (tt<alpha);
%             pdiff(:,i)=(tt<alpha);
%            end
     end

     warning ('off','all');
     try
         CC=bwareafilt(pdiff==1,1);
%          CC=ValidateRow(CC);
%          CC=ValidateRow(CC');
%          CC=bwareafilt(CC,1);
     catch
         CC=zeros(m,n);
     end
     
     mm=sum(sum(CC))/(m-1)/(n-1)
     thres=max([prc(j)*0.5/100;1/n]);
     if mm>thres
         %         sum(sum(CC))/(m-1)/(n-1)/prc(j)*100
         indAll=(CC==1);
         p=alpha;
         break;
     end

    
%      pp=zeros(n,1);
% %      sind=zeros(n,2);
%      for i=2:n
%          is=max(2,i-nn);
%          ie=min(n,i+nn);
%          tmp=sum(pdiff(:,is:ie),2);
%          tmp=find(tmp<ie-is+1);
%          tmp=[1;tmp;m+1];
%          tmp2=diff(tmp);
%          %ind=find(tmp2==max(tmp2),1,'last');
%          jm=max(tmp2);
% %          j=find(tmp2==jm,1,'last');
%          
%          %     if sum(P(tmp(j+1)-1,i-1:i+1)<alpha)==3
%          pp(i)=(jm-1)/(m-2);
% %          sind(i,:)=[tmp(j)+1, tmp(j+1)-1];
%          %     end
%          %      else
%          %           sind(i,:)=[1,1];
%          %      end
%      end
     
%      figure
%      plot(pp)
% max(pp)
     %thres=max(0.45,min(pp(n), mean(pp(2:n)))+3*std(pp(2:n)));
     %thres=max(0.3,min(mean(pp(n))+3*std(pp(2:n)), mean(pp(2:n))+3*std(pp(2:n))));
%      figure 
%      plot(pp)
%      thres=max(0.2,min(pp(n), mean(pp(2:n)))+3*std(pp(2:n)));
%      if max(pp)>=thres
%          p=alpha;
%          break;
%      end
end
% max(pp)
%

% [mm,i]=max(pp);
% % i
% % figure
% % plot(pp);
% %thres=0.4;
% thres=min([0.5, mean(pp(2:n))+4*std(pp(2:n))]);%,max(pp)]);
% % mm
% % thres
% p=1;
% if mm>thres
% %     mm
% %     thres
%     i=find(pp>thres);
%     p=1;
%     for t=1:length(i)
%         j=i(t);
%          pt=min(min(P(sind(j,1):sind(j,2),j)))/pp(j);
%          p=min(p,pt);
%     end
% % else
% %     if mm>0.2
% %     %p=min(P(end));%,median(P(P<1)));
% %     %j=n;
% %     %p=min(P(sind(j,1):sind(j,2),j))/pp(j);
% %     end
% end
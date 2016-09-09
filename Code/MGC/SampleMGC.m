function test=SampleMGC(P)
% Find the smooth regions in the p-value map, by considering the
% largest monotonically decreasing or increasing scales along
% the row or column p-values, but allowing small p-value increase or
% decrease as specified by tau.
[m,n]=size(P);
P2=P(2:end,2:end);
P2=P2(P2<0);
thres=1/min(m,n);
% st=std(P2);
st=max(norm(P2,'fro')/sqrt(length(P2)),0.005);
t1=prctile(P(P>0),90);
t2=7*st;
t3=max(P(end),0)+st;
t=max([t1,t2,t3]);
R=(P>t);
if mean(mean(R))>0
    warning('off','all');
    R=bwareafilt(R,1);
end

if mean(mean(R))<thres
    test=P(end);
else
    [k,l]=find(P==max(max(P(R))));
    ln=ceil(0.1*m);
    lm=ceil(0.1*n);
    test=0;
    for i=1:length(k)
        ki=k(i);
        li=l(i);
        if R(ki,li)==1
            left=max(2,li-ln);
            right=min(n,li+ln);
            upper=max(2,ki-lm);
            down=min(m,ki+lm);
            tmp1=min(P(upper:down,li));
            tmp2=min(P(ki,left:right));
            tmp=max(tmp1,tmp2);
            %tmp=P(upper:down,left:right);
           % tmp=median(tmp(tmp<2));
            if tmp>test
                test=tmp;
            end
        end
    end
%     if test<t
%         test=P(end);
%     end
end
% 
% %st=max(-min(P2),0.005);
% 
% %st=max(3*norm(P2(2:end,2:end),'fro')/sqrt((m-1)*(n-1)),0.03)
% % (prctile(P2(P2>0),95))/st
% % (prctile(P2(P2>0),90))/st
% % (prctile(P2(P2>-1),80))/st
% warning('off','all');
% R=(P>-mean(P2)+7*st);
% % R=(P>(P(end)+st));
% if mean(mean(R))>=thres
% %     tic
%     R=bwareafilt(R,1);
% %     [R,s]=Smooth(P,7*st);
%     %s
% %     toc
% end
% % % if mean(mean(tmp))>2/min(m,n)
% %     warning('off','all');
% %     tmp=bwareafilt(tmp,1);
% % end
% if mean(mean(R))<thres
%     %     tmp=bwareafilt(tmp,1);
%     test=P(end);
% else
%     test=max(P(R));
% %         test=median(P(R));
%     %     test=min(min(P(tmp)));
%     %     test=median((P2(tmp)));
%     % t1=max(median(P2,1));
%     % t2=max(median(P2,2));
%     % test=max([t1,t2]);
%     % test=prctile(P(P>0),90);
%     % %     cc=[];
% %     lm=ceil(0.02*m);
% %     ln=ceil(0.02*n);
% % % %         lm=1;
% % % %         ln=1;
% %     [k,l]=find(P==max(max(P)));
% % %     test=max(median(P(k,:)),median(P(:,l)));
% % %    [k,l]=find(P>=prctile(P(R),99));
% % %         test=max(max(P));
% %     test=0;
% %     for i=1:length(k)
% %         ki=k(i);
% %         li=l(i);
% %         if R(ki,li)==1
% %             left=max(2,li-ln);
% %             right=min(n,li+ln);
% %             upper=max(2,ki-lm);
% %             down=min(m,ki+lm);
% %             
% %             tmp=P(upper:down,left:right);
% %              tmp=median(tmp(tmp<2));
% % %             tmp1=min(P(ki,left:right));
% % %             tmp2=min(P(upper:down,li));
% % %             tmp=max(tmp1,tmp2);
% % %             tmp=max(tmp1,tmp2);
% %             if tmp>test
% %                 test=tmp;
% %             end
% %         end
% %     end
%     if test<P(end)+4*st
%         test=P(end);
%     end
% end

% function [R,s]=Smooth(P,thres)
% [m,n]=size(P);
% R2=cell(4,1);
% PD1=zeros(m,n);
% PD2=zeros(m,n);
% for i=2:m
%     PD1(i,2:end)=diff(P(i,:));
% end
% for i=2:n
%     PD2(2:end,i)=diff(P(:,i));
% end
% R=(P>thres);
% 
% R2{1}=(PD1>=0)&R;
% % R2{2}=(PD1<=0)&R;
% R2{2}=(PD2>=0)&R;
% % R2{4}=(PD2<=0)&R;
% R=zeros(m,n);
% s=0;
% for i=1:2
%     if mean(mean(R2{i}))>0
%     R2{i}=bwareafilt(R2{i},1);
%     t=mean(mean(R2{i}));
%     if t>s
%         R=R2{i};
%         s=t;
%     end
%     end
% end
% % figure
% % imagesc(R)
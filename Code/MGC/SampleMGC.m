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
st=max(norm(P2,'fro')/sqrt(length(P2)),0.004);
% t1=prctile(P(P>0),90);
t=7*st;
% t3=max(P(end),0)+st;
% t=max([t2,t3]);
% t=max(t2,t3);
R=(P>t);
% if mean(mean(R))>0
%     warning('off','all');
%     R=bwareafilt(R,1);
% end
    test=P(end);
    
if mean(mean(R))<2*thres
    test=P(end);
else
%     R=Smooth(P)&R;
%     if mean(mean(R))<=thres
%         test=P(end)
%     else
        [~, ~, ~, R] = FindLargestRectangles(Smooth(P,R), [0,0,1],[2,2]);
%             if mean(mean(R))<=thres
%         test=P(end)
%     else
        %     R=R2&R;
        if sum(sum(R))>0
            test=max(max(P(R)));
        end
        %         [k,l]=find(P>=max(max(P(R))));
        % %         [k,l]=find(R==1);
%         ln=ceil(0.1*m);
%         lm=ceil(0.1*n);
%         test=0;
%         %     test=max(max(P));
%         for i=1:length(k)
%             ki=k(i);
%             li=l(i);
%             if R(ki,li)==1
%                 left=max(2,li-ln);
%                 right=min(n,li+ln);
%                 upper=max(2,ki-lm);
%                 down=min(m,ki+lm);
%                             tmp1=median(P(upper:down,li));
%                             tmp2=median(P(ki,left:right));
%                             tmp=max(tmp1,tmp2);
% %                 tmp=P(upper:down,left:right);
%                 tmp=median(tmp(tmp<2));
%                 if tmp>test
%                     test=tmp;
%                 end
%             end
%         end
% %                 test=median(P(R));
%             if test<t
%                 test=P(end)
%             end
%     end
end

function [R]=Smooth(P,R)
if nargin<2
    R=ones(size(P));
end
[m,n]=size(P);
R2=cell(4,1);
PD1=zeros(m,n);
PD2=zeros(m,n);
for i=2:m
    PD1(i,2:end)=diff(P(i,:));
end
for i=2:n
    PD2(2:end,i)=diff(P(:,i));
end

R2{1}=(PD1>=0)&R;
R2{2}=(PD1<=0)&R;
R2{3}=(PD2>=0)&R;
R2{4}=(PD2<=0)&R;
R=false(m,n);
for i=1:4
    if sum(sum(R2{i}))>0
        R2{i}=bwareafilt(R2{i},1);
        R= (R | R2{i});
%         t=mean(mean(R2{i}));
%         if t>s
%             R=R2{i};
%             s=t;
%         end
    end
end
% % figure
% % imagesc(R)
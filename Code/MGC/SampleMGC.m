function test=SampleMGC(P)
% Find the smooth regions in the p-value map, by considering the
% largest monotonically decreasing or increasing scales along
% the row or column p-values, but allowing small p-value increase or
% decrease as specified by tau.
[m,n]=size(P);
P2=P(2:end,2:end);
P2=P2(P2<0);
thres=1/min(m,n);
t=3.5*max(norm(P2,'fro')/sqrt(length(P2)),0.01);
R=(P>t);
test=P(end);

if mean(mean(R))>=min(2*thres,0.05)
    R=Monotone(P,R);
    %     [~, ~, ~, R] = FindLargestRectangles(Smooth(P,R), [0,0,1],[2,2]);
    if mean(mean(R))>0
        [k,l]=find((P>=max(max(P(R))))&(R==1));
        ln=ceil(0.1*m);
        lm=ceil(0.1*n);
        for i=1:length(k)
            ki=k(i);
            li=l(i);
            left=max(2,li-ln);
            right=min(n,li+ln);
            upper=max(2,ki-lm);
            down=min(m,ki+lm);
            tmp1=min(P(upper:down,li));
            tmp2=min(P(ki,left:right));
            tmp=max(tmp1,tmp2);
            if tmp>test
                test=tmp;
            end
        end
    end
end

function [R]=Monotone(P,R)
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
warning('off','all');
for i=1:4
    t=sum(sum(R2{i}));
    if t>0
        R2{i}=bwareafilt(R2{i},1);
        R= (R | R2{i});
    end
end
function test=SampleMGC(P)
% Find the smooth regions in the p-value map, by considering the
% largest monotonically decreasing or increasing scales along
% the row or column p-values, but allowing small p-value increase or
% decrease as specified by tau. 
[m,n]=size(P);
P2=P(2:end,2:end);
% st=max(norm(P2(P2<0),'fro')/sqrt(length(P2(P2<0))),0.03);
st=max(norm(P2(2:end,2:end),'fro')/sqrt((m-1)*(n-1)),0.03);
% (prctile(P2(P2>0),95))/st
% (prctile(P2(P2>0),90))/st
% (prctile(P2(P2>-1),80))/st
if (prctile(P2(P2>-1),80))/st<1
    test=P(end);
else
%     cc=[];
    lm=ceil(0.01*m);
    ln=ceil(0.01*n);
%      [k,l]=find(P==max(max(P)));
    [k,l]=find(P>=prctile(P(P<2),99));
    test=max(max(P));
%     test=0;
    for i=1:length(k)      
        ki=k(i);
        li=l(i);

        left=max(2,li-ln);
        right=min(n,li+ln);
        upper=max(2,ki-lm);
        down=min(m,ki+lm);
        tmp=P(upper:down,left:right);
        tmp=median(tmp(tmp<1));
        if tmp<test
            test=tmp;
        end
% %         RT(upper:down,left:right)=0;
%         tmp=tmp(tmp>-1);
%         if size(tmp,2)~=1;
%             tmp=tmp';
%         end
%         cc=[cc;tmp];
    end
end
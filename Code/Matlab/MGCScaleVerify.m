function [ind]=MGCScaleVerify(V)
% An auxiliary function to verify and estimate the MGC optimal scale

VN=V(2:end,2:end);
k=Verify(VN)+1;
l=Verify(VN')+1;
ind=find(V==min(min(V(k,l))),1,'last');

function k=Verify(VN)
thres=0.05;
[m,~]=size(VN);
lim=2;

k=m;
if median(VN)<thres
    k=1:m;
end
% 
% rowTmp=median(VN,2)';
% [~,indK]=sort(rowTmp,'ascend');
% 
% 
% for l=1:lim
%     switch l
%         case 1
%             tmpThres=thres;
% %         case 2
% %             tmpThres=thres/2;
%         case 2 
%             tmpThres=thres/5;
%     end
% %     tmpThres=thres/lim*(lim-l+1);
%     
%     cE=find(rowTmp>tmpThres);
%     cS=[0 cE];
%     cE=[cE m+1];
%     c=max(cE-cS)-1;
%     if c>=m*tmpThres/thres/2
%         tt=find((cE-cS)==max(cE-cS),1,'last');
%         k=cS(tt)+1:cE(tt)-1;
%         break;
%     end
%     
% %     len=indK(1:ceil(m*tmpThres/thres));
% %     if median(rowTmp(len))<=tmpThres
% %         k=len;
% %         break;
% %     end
% end
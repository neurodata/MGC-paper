function  [p,pAll,test,testAll,indAll]=MGCPermutationTest(A,B,rep,option)
% Author: Cencheng Shen
% This function tests independent between two data sets, using MGC by a random permutation test.
%
% The inputs are:
% two n*n distance matrices A & B,
% a parameter rep to specify the number of random permutations,
% an option to specify which global test to use, set to 1,2,3 for dcorr / mcorr / Mantel.
%
% The outputs are:
% the estimated MGC p-value, the p-values of all local tests, 
% the estimated MGC test statistic, all local test statistics, and the estimated optimal scale. 

if nargin<3
    rep=1000; % use 1000 random permutations by default
end
if nargin<4
    option=2;  % use mcorr by default
end
[m,n]=size(A);

% calculate all local correlations between the two data sets
testAll=LocalCorr(A,B,option);
 test=SmoothRegions(testAll);
% calculate the local correlations under permutation, to yield the p-values of all observed local correlations
for r=1:rep
    % use random permutations on the second data set
    per=randperm(n);
    BN=B(per,per);
    tmp=LocalCorr(A,BN,option);
    tmp2=SmoothRegions(tmp);
    if r==1
        pAll=(tmp>=testAll)/rep;
        p=(tmp2>=test)/rep;
    else
        pAll=pAll+(tmp>=testAll)/rep;
        p=p+(tmp2>=test)/rep;
    end
end
if pAll(end)<0.025
    p=pAll(end);
    test=testAll(end);
end
% % set the p-values of local corr at rank 0 to maximum, since they should not be used
% if (sum(sum(pAll<1/rep))>0)
%     pAll=pAll+1/rep;
% end
% pAll(pAll>1)=1;
% pAll(1,:)=1;pAll(:,1)=1;
warning('off','all');
[~,~,~,indAll]=FindLargestRectangles((pAll<=p), [0 0 1],[2,2]);
indAll=find(indAll==1);
[m,n]=size(pAll);
if pAll(end)<=p && sum(indAll==m*n)==0
    indAll=[indAll;m*n];
end

function test=SmoothRegions(P)
% Find the smooth regions in the p-value map, by considering the
% largest monotonically decreasing or increasing scales along
% the row or column p-values, but allowing small p-value increase or
% decrease as specified by tau. 
[m,n]=size(P);
[k,l]=find(P>=prctile(P(P<2),95));
% [k,l]=find(P==max(max(P)));
cc=[];
lm=ceil(0.025*m);
ln=ceil(0.025*n);
RT=ones(m,n);
for i=1:length(k)
    ki=k(i);
    li=l(i);
%     tmp=max([median(P(ki,:));median(P(:,li))]);
    left=max(2,li-ln);
    right=min(n,li+ln);
    upper=max(2,ki-lm);
    down=min(m,ki+lm);
    tmp=P(upper:down,left:right);
    RT(upper:down,left:right)=0;
    tmp=tmp(tmp>-1);
    if size(tmp,2)~=1;
        tmp=tmp';
    end
   cc=[cc;tmp];
end
[h,edges]=histcounts(cc,min(ceil(length(cc)/10),10)); % partition the smooth p-values into bins
[~,ii]=max(h); % find the bin with most counts
test=edges(ii); % use the largest p-value in that bin for MGC
% R=bwareafilt(P>=test,1);
% [~,~,~,R]=FindLargestRectangles((P>=test), [0 0 1],[2,2]);
% mean(mean(R))
% P2=P;
% P2(RT==1)=0;
% R=FindLargestSquares((P>=test));
%  max(max(R))^2/m/n
% max(max(R))^2/m/n
% if max(max(R))^2/m/n<0.008 || test<P(end)
%     test=P(end);
% end

% 
% t1=zeros(m,1);
% t2=zeros(n,1);
% 
% PD=cell(2,1);
% PD{1}=zeros(m,n); % store the p-value changes within rows
% PD{2}=zeros(m,n); % store the p-value changes within columns
% for i=2:m
%     tt=P(i,:);
%     PD{1}(i,2:end)=diff(tt);   
% end
% for i=2:n
%     tt=P(:,i);
%     PD{2}(2:end,i)=diff(tt);
% end
% 
% RC=cell(4,1);
% 
% RC{1}=(PD{1}<=tau); % check monotone decreasing , but also allows small p-value increase no more than tau
% RC{2}=(PD{1}>=-tau); % check monotone increasing , but also allows small p-value decrease no more than tau
% RC{3}=(PD{2}<=tau); % repeat for column changes
% RC{4}=(PD{2}>=-tau);
% 
% % kk=ceil(0.02*m);
% % ll=cel(0.02*n);
% for i=2:m
% %     b1=max(2,i-kk);
% % %     b2=min(m,i+kk);
% %     tmp=find(sum(RC{2}(b1:b2,2:end),2)==(2*kk+1));
%     tmp=find(RC{2}(i,2:end)==0);
%     tmp=diff(tmp);
%     t1(i)=max(tmp)/m;
% end
% 
% for j=2:n
% %      b1=max(2,i-kk);
% %     b2=min(n,i+ll);
%     tmp=find(RC{4}(2:end,j)==0);
%     tmp=diff(tmp);
%     t2(j)=max(tmp)/n;
% end
% 
% k=find(t1>prctile(t1,90));
% l=find(t2>prctile(t2,90));
% test=P(k,l);
% test=test(test>-1);
% if max(t1)<0.1 || max(t2)<0.1
%     test=P(end);
% else
%     test=median(test);
% end


% % tau=0;
% sumC=median(P(2:end,2:end),1);
% sumR=median(P(2:end,2:end),2);
% 
% % test=prctile(sumC,95);
% % l=find(sumC>test)+1;
% % test=prctile(sumR,95);
% % k=find(sumR>test)+1;
% gind=0;
% 
% % if sum(k==m)>0 && sum(l==n)>0
% %     gind=1;
% % end
% % l=find(sumC==max(sumC))+1;
% % k=find(sumR==max(sumR))+1;
% % if mean(sumC(end)<sumC)<0.05 && mean(sumR(end)<sumR)<0.05
% %     gind=1;
% % else
% %     gind=0;
% % end
% 
% % 
% % % test=max(max(P(k,l)));
% % tmp=P(k,l);
% % tmp=tmp(tmp>-1);
% % if length(tmp)>1
% %     %     tmp=P(k,l);
% %     %     tmp=tmp(tmp>-1);
% %     [h,edges]=histcounts(tmp,10); % partition the smooth p-values into bins
% %     [~,ii]=max(h); % find the bin with most counts
% %     test=edges(ii); % use the largest p-value in that bin for MGC
% % else
% %     test=max(tmp);
% % end
% test=P(P>-1);
% test=prctile(test,95);
% % [k,l]=find(P>=test);
% % 
% % tmp=P(k,2:end);
% % tmp=tmp(tmp>-1);
% % tmp2=P(2:end,l)';
% % tmp2=tmp2(tmp2>-1);
% % cc=[tmp;tmp2];
% % cc=P(k,l);
% % 
% % lm=ceil(0.02*m);
% % ln=ceil(0.02*n);
% % % test=max(max(P));% % test=0;
% % cc=[];
% % for i=1:length(k)
% %     ki=k(i);
% %     li=l(i);
% % %     tmp=max([median(P(ki,:));median(P(:,li))]);
% %     left=max(2,li-ln);
% %     right=min(n,li+ln);
% %     upper=max(2,ki-lm);
% %     down=min(m,ki+lm);
% %     tmp=P(upper:down,left:right);
% %     tmp=tmp(tmp>-1);
% %     if size(tmp,2)~=1;
% %         tmp=tmp';
% %     end
% %    cc=[cc;tmp];
% % %     tmp=median(tmp);
% % %     tmp=min(tmp);
% % %     [h,edges]=histcounts(tmp,10); % partition the smooth p-values into bins
% % %     [~,ii]=max(h); % find the bin with most counts
% % %     tmp=edges(ii); % use the largest p-value in that bin for MGC
% % 
% % %     tmp=min(min(P(upper:down,left:right)));
% % %     tmp=mean(mean(P(upper:down,left:right)));
% % %     if tmp>test
% % %         test=tmp;
% % %     end
% % end
% [h,edges]=histcounts(cc,10); % partition the smooth p-values into bins
% [~,ii]=max(h); % find the bin with most counts
% test=edges(ii); % use the largest p-value in that bin for MGC
% 
% % if mean(mean(P(end)<P))<0.05
% %     test=P(end);
% % end


% % % 
% tau=0;
% PD=cell(2,1);
% PD{1}=zeros(m,n); % store the p-value changes within rows
% PD{2}=zeros(m,n); % store the p-value changes within columns
% for i=2:m
%     tt=P(i,:);
%     PD{1}(i,2:end)=diff(tt);
% end
% for i=2:n
%     tt=P(:,i);
%     PD{2}(2:end,i)=diff(tt);
% end
% 
% RC=cell(4,1);
% 
% RC{1}=(PD{1}<=tau); % check monotone decreasing , but also allows small p-value increase no more than tau
% RC{2}=(PD{1}>=-tau); % check monotone increasing , but also allows small p-value decrease no more than tau
% RC{3}=(PD{2}<=tau); % repeat for column changes
% RC{4}=(PD{2}>=-tau); 
% 
% % CC=[mean(P,1)';mean(P,2)];
% % % C2=median(P,2);
% % test=max(CC);
% R=zeros(m,n);
% sz=0;
% for i=1:4
%     if i<3
%         [R2]=verify(RC{i});
%     else
%         [R2]=verify(RC{i}');
%         R2=R2';
%     end
%     t=mean(mean(R2));
%     if t>sz
%         sz=t;
%         R=R2;
%     end
% end
% figure
% imagesc(R)
% R=bwareafilt(R==1,1);
% % tic
% % for i=1:4
% %     RC{i} = FindLargestSquares(RC{i}); % find the largest rectangle within each (approximately) monotonically decreasing / increasing region
% %     R=bsxfun(@max,R,RC{i});
% % end
% 
% %R=(RC{1} | RC{2} | RC{3} | RC{4}); % combine all four rectangles together into the smooth regions
% 
% %[~,~,~,R2]=FindLargestRectangles(RC{2}, [0 0 1],[2,2]);
% % tmp=FindLargestSquares(RC{2});
% % sz=max(max(R))
% mean(mean(R))
% if sz==0
%     test=P(end);
% else
% %     [k,l]=find(R==sz);
% %     R2=zeros(m,n);
% %     for i=1:length(k)
% %         R2(k(i):k(i)+sz-1,l(i):l(i)+sz-1)=1;
% %     end
%     % test=max(P(R2));
%     %
%     [h,edges]=histcounts(P(R==1),min(ceil(sum(sum(R))/10),100)); % partition the smooth p-values into bins
%     [~,ii]=max(h); % find the bin with most counts
%     test=edges(ii); % use the largest p-value in that bin for MGC
% end
% 
% function [R2,s]=verify(R)
% [m,n]=size(R);
% R2=zeros(m,n);
% kl=ceil(0.025*m);
% s=0;
% i=2;
% while (i<=m)
%     k1=max(2,i-kl);
%     k2=min(m,i+kl);
%     tmp=find(sum(R(k1:k2,:),1)<(2*kl+1));
%     [tmpM,ind]=max(diff(tmp));
%     if tmpM/n>0.05
%         ind=tmp(ind)+1:tmp(ind+1)-1;
%         R2(k1:k2,ind)=1;
%         i=k2+1;
%         s=(tmpM-1)/n;
%     else
%         i=i+1;
%     end
% end

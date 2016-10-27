function plot_panel1(F,x,y,R2)

siz=size(x);
dim=siz(2);
n=siz(1);
if dim>1
    xnew=zeros(n,1);ynew=zeros(n,1);
    for i=1:n
        %if group(i)==0
        tmp=find(R2(:,i)==1)';
        %         tmp=[tmp];
        %end
        if length(tmp)>1
            [~,~,~,xt,yt]=canoncorr(x(tmp,:),y(tmp,:));
            tmp2=find(tmp==i);
            xnew(tmp2)=xt(tmp2,1);ynew(i)=yt(tmp2,1);
        end
        %         plot(x(tmp),yest(tmp),'-','Color',loca,'linewidth',3);
    end
    x=xnew;y=ynew;
end
% ax=subplot(s,t,1);
subplot('Position',F.pos);
hold all
set(groot,'defaultAxesColorOrder',F.map2);
plot(x,y,'.','MarkerSize',F.mkSize,'Color',F.gray);
% if F.type==1 && noise==1
%     xlabel('Cloud Shape','FontSize',F.fontSize,...
%         'Units', 'normalized','Position', [-0.01, -0.2], 'HorizontalAlignment', 'left')
%     ylabel('Ground Wetness','FontSize',fontSize, ...
%         'Units', 'normalized', 'Position', [-0.06 0], 'HorizontalAlignment', 'left')
% else
xlabel('$x$','FontSize',F.fontSize2+6,'Interpreter','latex',...
    'Units', 'normalized','Position', [-0.01, -0.2], 'HorizontalAlignment', 'left')
ylabel('$y$','FontSize',F.fontSize2+6,'Interpreter','latex', ...
    'Units', 'normalized', 'Position', [-0.06 0], 'HorizontalAlignment', 'left')
% end

% if F.type>5
%     [I,J]=ind2sub([n,n],find(C_MGC>0.1,1,'first'));
%     J2=find(mcorrH(J,:)<0,1,'last');
% else
% end

for ind=[1,2,3]; %length(F.id)
    text(x(F.id(ind))+F.hs,y(F.id(ind))+F.hy(ind),num2str(ind),'fontsize',F.fontSize,'color',F.col)
    plot(x(F.id(ind)),y(F.id(ind)),'.','MarkerSize',F.mkSize,'Color',F.col);
end

% tname=CorrSimuTitle(F.type);
% findex=strfind(tname,'.');
% tname=tname(findex+1:end);
% xlim([min(x)-0.2, max(x)]);
% ylim([min(y)-0.2, max(y)]);

set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
axis('square')

if F.tit
    tit1=strcat('0', F.AB ,'. Sample Data');
    title([{tit1}; {' '}], 'Units', 'normalized', ...
        'Position', [0 1.1], 'HorizontalAlignment', 'left','FontSize',F.tfs);
end
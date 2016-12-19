function plot_panel2(F,C,D)

siz=size(C);
n=siz(1);

C2=reshape(C,n^2,1);
D2=reshape(D,n^2,1);

ax=subplot('Position',F.pos2);
% ax=subplot(s,t,2);
RC=DistRanks(C);
RD=DistRanks(D)';
RC=(RC<=F.Xmax+1);
RD=(RD<=F.Ymax+1);
ind1=reshape(RC&RD,n^2,1);
hold on
set(groot,'defaultAxesColorOrder',F.map2);

if F.type==1
I=4;
else
    I=21;
end
%%% Here I subsample one point
C3=reshape(C(I,:),n,1);
D3=reshape(D(I,:),n,1);
ind2=reshape(RC(I,:)&RD(:,I)',n,1);

plot(C3,D3,'.','MarkerSize',6,'Color',F.gray);
plot(C3(ind2==1),D3(ind2==1),'o','MarkerSize',4,'Color',F.loca);
% This plots all points
% plot(C2,D2,'.','MarkerSize',6,'Color',F.gray);
% plot(C2(ind1==1),D2(ind1==1),'o','MarkerSize',4,'Color',F.loca);



%%% useless codes
% tmpX=C2(ind1==0);
% tmpY=D2(ind1==0);
% t= 0:pi/10:2*pi;
% for i=1:length(tmpX)
% %     i
% pb=patch((0.01*sin(t)+ tmpX(i)),(0.01*cos(t)+tmpY(i)),F.gray,'edgecolor','none');
% alpha(pb,0.1);
% end
% plot(C2(ind1==1),D2(ind1==1),'+','MarkerSize',4,'Color',F.loca);

x12=sub2ind([n,n], F.id(1),F.id(2));
x23=sub2ind([n,n], F.id(2),F.id(3));
text(C2(x12)+0.02,D2(x12),'(1, 2)','fontsize',F.fontSize,'color',F.col)
plot(C2(x12),D2(x12),'.','MarkerSize',8,'Color',F.col);

text(C2(x23)+0.02,D2(x23),'(2, 3)','fontsize',F.fontSize,'color',F.col)
plot(C2(x23),D2(x23),'.','MarkerSize',8,'Color',F.col);
hold off
alpha(0.1)

xlim([min(min(C)),1]);
ylim([min(min(D)),1]);
warning('off','all')
if F.type==1
    xlabel('$$d_{x}(x_i,x_j)$$','FontSize',F.fontSize2+2,...
        'Units', 'normalized','Position', [0.5, -0.02], 'HorizontalAlignment', 'center','Interpreter','latex');
    ylabel('$$d_{y}(y_i,y_j)$$','FontSize',F.fontSize2+2, ...
        'Units', 'normalized', 'Position', [-0.28 0.5], 'HorizontalAlignment', 'center','Interpreter','latex');
    set(gca,'XTick',[0,1],'YTick',[0,1],'FontSize',F.fontSize); % Remove x axis tick
else
    set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
end
txt2 = strcat('\color[rgb]{0.5 0.5 0.5} c(Dcorr)=', num2str(round(100*F.tA(end))/100));
title({txt2},'FontSize',F.tfs,'interpreter','tex');

% if F.tit
%     tit1=strcat('1', F.AB ,'. Pairwise Distances');
%     title([{tit1}; {' '}], 'Units', 'normalized', ...
%         'Position', [0 1.1], 'HorizontalAlignment', 'left','FontSize',F.tfs);
% end
% set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
set(gca,'FontSize',F.fontSize); % Remove x axis tick
axis('square')
pos2 = get(ax,'position');
pos2(3:4) = F.pos(3:4);
set(ax,'position',pos2);

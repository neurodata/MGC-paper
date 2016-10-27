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
plot(C2(ind1==0),D2(ind1==0),'.','MarkerSize',6,'Color',F.gray);
plot(C2(ind1==1),D2(ind1==1),'+','MarkerSize',4,'Color',F.loca);

x12=sub2ind([n,n], F.id(1),F.id(2));
x23=sub2ind([n,n], F.id(2),F.id(3));
text(C2(x12)+F.hs,D2(x12)+F.hy(1),'(1, 2)','fontsize',F.fontSize,'color',F.col)
plot(C2(x12),D2(x12),'.','MarkerSize',F.mkSize/2,'Color',F.col);

text(C2(x23)+F.hs,D2(x23)+F.hy(1),'(2, 3)','fontsize',F.fontSize,'color',F.col)
plot(C2(x23),D2(x23),'.','MarkerSize',F.mkSize/2,'Color',F.col);
hold off

% tname=CorrSimuTitle(F.type);
% findex=strfind(tname,'.');
% tname=tname(findex+1:end);
xlim([min(C2), max(C2)]);
ylim([min(D2), max(D2)]);
warning('off','all')
xlabel('Distance$$_{x}(x_i,x_j)$$','FontSize',F.fontSize2+4,...
    'Units', 'normalized','Position', [-0.01, -0.2], 'HorizontalAlignment', 'left','Interpreter','latex');
ylabel('Distance$$_{y}(y_i,y_j)$$','FontSize',F.fontSize2+4, ...
    'Units', 'normalized', 'Position', [-0.06 0], 'HorizontalAlignment', 'left','Interpreter','latex');

if F.tit
    tit1=strcat('1', F.AB ,'. Pairwise Distances');
    title([{tit1}; {' '}], 'Units', 'normalized', ...
        'Position', [0 1.1], 'HorizontalAlignment', 'left','FontSize',F.tfs);
end
set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
%pos=[nan, nan, width, height];
axis('square')
pos2 = get(ax,'position');
pos2(3:4) = F.pos(3:4);
set(ax,'position',pos2);

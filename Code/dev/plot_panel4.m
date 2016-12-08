function plot_panel4(F,tN,tA,k,l,testN,mcorrH,A,B,test,pMLocal,pMGC)

n=size(tA,2);

% Null Distributions
ax=subplot('Position',F.pos4);
minp=min([min(tN(:,n,n)),min(tN(:,k,l)),tA(k,l),tN(n,n)]);
minp=floor(minp*10)/10;
maxp=max([max(tN(:,n,n)),max(tN(:,k,l)),tA(k,l),tN(n,n)]);
maxp=ceil(maxp*10)/10;
p=tN(:,k,l);
[~,~]=ksdensity(p,'support',[-1,1]);
p=tN(:,n,n);
[f,xi]=ksdensity(p,'support',[-1,1]);
p=testN;
[f2,xi2]=ksdensity(p,'support',[-1,1]);

hold on
plot(xi,f,'.-','LineWidth',4,'Color',F.glob);
plot(xi2,f2,'.-','LineWidth',4,'Color',F.mgc);
set(gca,'FontSize',F.fontSize);
x1=sum(sum(mcorrH))/norm(A,'fro')/norm(B,'fro');
x1=round(x1*100)/100;
% x2=sum(sum(C_MGC))/norm((A_MGC-mean(mean(A_MGC))),'fro')/norm((B_MGC--mean(mean(B_MGC))),'fro');
% x2=round(x2*100)/100;
x3=round(test*100)/100;
plot(x1,0.1,'*','MarkerSize',12,'Color',F.glob,'linewidth',2);
plot(x3,0.1,'*','MarkerSize',12,'Color',F.mgc,'linewidth',2);

% ind=find(xi>x1,1,'first');
y1=max(f)+2;
% y2 = max(f1)+2;
y3 = 5;
txt1 = strcat('$$p(c) =', num2str(pMLocal(end)),'$$');
txt3 = strcat('$$p(\hat{c}^{*}) = ', num2str(pMGC),'$$');
c=text(x3,y3,txt3,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',F.mgc,'Interpreter','latex');
set(c,'FontSize',F.fontSize);
set(gca,'XTick',x3+0.1,'TickLength',[0 0],'XTickLabel',x3);
a=text(x1,y1,txt1,'VerticalAlignment','bottom','HorizontalAlignment','left','Color',F.glob,'Interpreter','latex');
set(a,'FontSize',F.fontSize);
if abs(x1-x3)>0.02
    set(gca,'XTick',sort([x3+0.05,x1+0.05]),'TickLength',[0 0],'XTickLabel',sort([x3,x1]));
end
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',F.fontSize);
ylim([0 y1+10]);

xlim([minp,min(maxp+0.1,1)]);
if F.type<5
    xlim([minp,min(maxp+0.4,1)]);
end
xlabel('Test Statistic','FontSize',F.fontSize2,...
    'Units', 'normalized','Position', [-0.01, -0.20], 'HorizontalAlignment', 'left');
ylabel('Density','FontSize',F.fontSize2, ...
    'Units', 'normalized', 'Position', [-0.10, 0], 'HorizontalAlignment', 'left');
set(gca,'YTick',[])

if F.tit
    tit1=strcat('3', F.AB ,'. Null Distributions');
    title([{tit1}; {'& P-Values'}],'FontSize',F.tfs, ...
        'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left')
end
axis('square')
hold off
pos2 = get(ax,'position');
pos2(3:4) = F.pos(3:4);
set(ax,'position',pos2);
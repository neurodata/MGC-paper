function plot_panel5(F,pMLocal)

n=size(pMLocal,2);

ax=subplot('Position',F.pos5);
% ax=subplot(s,t,5);
hold on
set(groot,'defaultAxesColorOrder',F.map1);
kmin=2;
ph=pMLocal(kmin:n,kmin:n)';
imagesc(log(ph));

% draw boundary around optimal scale
% indP=optimalInd;
% [J,I]=ind2sub(size(pMLocal),indP);
% Ymin=min(I)-1;
% Ymax=max(I)-1;
% Xmin=min(J)-1;
% Xmax=max(J)-1;
lw=3;
plot([F.Xmin,F.Xmin],[F.Ymin,F.Ymax],'g','linewidth',lw)
plot([F.Xmax,F.Xmax],[F.Ymin,F.Ymax],'g','linewidth',lw)
plot([F.Xmin,F.Xmax],[F.Ymin,F.Ymin],'g','linewidth',lw)
plot([F.Xmin,F.Xmax],[F.Ymax,F.Ymax],'g','linewidth',lw)
xlim([2,n]);
ylim([2,n]);
%     imagesc(k,l,1);
hold off

set(gca,'FontSize',F.fontSize)
set(gca,'YDir','normal')
cmap=F.map4;
colormap(ax,flipud(cmap));
%ceil(max(max(ph))*10)/10
% caxis([0 1]);
cticks=[0.001, 0.01, 0.1, 0.5];
h=colorbar('Ticks',log(cticks),'TickLabels',cticks);%,'location','eastoutside');
set(h,'FontSize',F.fontSize);
axpos = ax.Position;
hpos=h.Position;
hpos(3)=0.5*hpos(3);
hpos(1)=hpos(1)+0.023;
h.Position=hpos;
ax.Position = axpos;
xlim([1 n-1]);
ylim([1 n-1]);
set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',F.fontSize);
xlabel('X Scales','FontSize',F.fontSize2,...
    'Units', 'normalized','Position', [-0.010, -0.20], 'HorizontalAlignment', 'left');
ylabel('Y Scales','FontSize',F.fontSize2, ...
    'Units', 'normalized', 'Position', [-0.22 -0.02], 'HorizontalAlignment', 'left');

if F.tit
    tit1=strcat('4', F.AB ,'. Multiscale P-Value');
    title([{tit1}; {'Map & Optimal Scales'}],'FontSize',F.tfs, ...
        'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left','color','g')
end

axis('square')
pos2 = get(ax,'position');
pos2(3:4) = F.pos(3:4);
set(ax,'position',pos2);

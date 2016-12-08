function plot_panel5(F,pMLocal,pMGC)

[m,n]=size(pMLocal);

ax=subplot('Position',F.pos4);
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
lw=2;
plot([F.Xmin,F.Xmin],[F.Ymin,F.Ymax],'g','linewidth',lw)
plot([F.Xmax,F.Xmax],[F.Ymin,F.Ymax],'g','linewidth',lw)
plot([F.Xmin,F.Xmax],[F.Ymin,F.Ymin],'g','linewidth',lw)
plot([F.Xmin,F.Xmax],[F.Ymax,F.Ymax],'g','linewidth',lw)
plot(F.k,F.l,'gs','markerSize',3,'MarkerFaceColor','g')
xlim([2,n]);
ylim([2,n]);
%     imagesc(k,l,1);

set(gca,'FontSize',F.fontSize)
set(gca,'YDir','normal')
cmap=F.map4;
colormap(ax,flipud(cmap));
%ceil(max(max(ph))*10)/10
cticks=[0.001, 0.01, 0.1, 0.5];
caxis(log([0.001 0.2]));
if F.type==8
    h=colorbar('Ticks',log(cticks),'TickLabels',cticks);%,'location','eastoutside');
    set(h,'FontSize',F.fontSize);
%     axpos = ax.Position;
%     hpos=h.Position;
%     % hpos(3)=0.5*hpos(3);
%     hpos(4)=0.5*hpos(4);
%     hpos(2)=hpos(2)-0.025;
%     % hpos(1)=hpos(1)+0.023;
%     h.Position=hpos;
%     ax.Position = axpos;
    set(h,'Box','off');
end
xlim([1 n-1]);
ylim([1 n-1]);

% % add MGC p-value
% txt1 = strcat('p(MGC) =', {' '},num2str(pMGC));
% if F.k<m/2
%     ph='left';
%     poh=1;
% else
%     ph='right';
%     poh=-1;
% end
% if F.l<n/2
%     pv='bottom';
%     pov=1;
% else
%     pv='top';
%     pov=-1;
% end
%
% a=text(F.k+poh,F.l+pov,txt1,'VerticalAlignment',pv,'HorizontalAlignment',ph,'Color','g');
% set(a,'FontSize',F.fontSize);
%
% % add global p-value
% txt1 = strcat('p(Dcorr) =', {' '}, num2str(pMLocal(m,n)));
plot(m-1,n-1,'s','markerSize',3,'MarkerFaceColor',F.glob,'Color',F.glob)
% a=text(m-1,n,txt1,'VerticalAlignment','bottom','HorizontalAlignment','right','Color',F.glob);
% set(a,'FontSize',F.fontSize);

hold off

if F.type==1
    xlabel('# X Neighbors','FontSize',F.fontSize2+2,...
        'Units', 'normalized','Position', [0.5, -0.2], 'HorizontalAlignment', 'center');
    ylabel('# Y Neighbors','FontSize',F.fontSize2+2, ...
        'Units', 'normalized', 'Position', [-0.28 0.5], 'HorizontalAlignment', 'center');
    set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',F.fontSize);
else
    set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
end

% if F.tit
%     tit1=strcat('3', F.AB ,'. Multiscale P-Value');
%     title([{tit1}; {'Map & Optimal Scales'}],'FontSize',F.tfs, ...
%         'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left','color','g')
% end
txt1 = strcat('\color[rgb]{0 1 0}p(MGC)=',num2str(pMGC));
txt2 = strcat('\color[rgb]{0.5 0.5 0.5}p(Dcorr)=', num2str(pMLocal(m,n)));
title({txt1,txt2},'FontSize',F.tfs-2,'interpreter','tex');

axis('square')
pos2 = get(ax,'position');
pos2(3:4) = F.pos(3:4);
set(ax,'position',pos2);

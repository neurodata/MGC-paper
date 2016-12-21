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
lw=1.5;
plot([F.Xmin,F.Xmin],[F.Ymin,F.Ymax],'g','linewidth',lw)
plot([F.Xmax,F.Xmax],[F.Ymin,F.Ymax],'g','linewidth',lw)
plot([F.Xmin,F.Xmax],[F.Ymin,F.Ymin],'g','linewidth',lw)
plot([F.Xmin,F.Xmax],[F.Ymax,F.Ymax],'g','linewidth',lw)
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
    set(h,'Box','off');
end
xlim([1 n]);
ylim([1 n]);

% plot scale points
plot(m,n,'.','markerSize',15,'MarkerFaceColor',F.glob,'Color',F.glob)
plot(F.k,F.l,'g','marker','o','markerSize',6,'linewidth',1)

hold off

if F.type==1
    xlabel('# X Neighbors','FontSize',F.fontSize2+2,...
        'Units', 'normalized','Position', [0, -0.2], 'HorizontalAlignment', 'left');
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
txt1 = strcat('\color[rgb]{0 1 0}p(MGC) = ',num2str(pMGC));
txt2 = strcat('\color[rgb]{0.5 0.5 0.5}p(Dcorr) = ', num2str(pMLocal(m,n)));
title({txt1,txt2},'FontSize',F.tfs,'interpreter','tex');

axis('square')
pos2 = get(ax,'position');
pos2(3:4) = F.pos(3:4);
set(ax,'position',pos2);

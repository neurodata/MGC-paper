function clean_panel(ax,map2,pos,id,n,col,fontSize)

colormap(ax,map2)
pos2 = get(ax,'position');
pos2(3:4) = [pos(3:4)];
set(ax,'position',pos2);
% plot([id(1)-0.5 id(4)-0.5],[id(1)-0.5 id(1)-0.5],'color',col)
% plot([id(1)-0.5 id(4)-0.5],[id(4)+0.5 id(4)+0.5],'color',col)
% plot([id(1)-0.5 id(1)-0.5],[id(1)-0.5 id(4)+0.5],'color',col)
% plot([id(4)-0.5 id(4)-0.5],[id(1)-0.5 id(4)+0.5],'color',col)
axis([1 n 1 n])
set(gca,'visible','off')
set(gca, 'visible', 'off')
set(gca,'FontSize',fontSize)
set(findall(gca, 'type', 'text'), 'visible', 'on')
set(gca,'XTick',[],'YTick',[]); % Remove y axis ticks
axis('square')
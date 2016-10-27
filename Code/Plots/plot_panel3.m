function plot_panel3(F,x,y,R2)
ax=subplot('Position',F.pos3);
hold on

n=size(x,1);

regressionLine=1; %% Regression line trial
if regressionLine==1
    yest=nan(n,1);
    for i=1:n
        %if group(i)==0
        tmp=find(R2(:,i)==1)';
        tmp=[tmp i];
        %end
        tmp2=[ones(length(tmp),1) x(tmp)];
        beta=tmp2 \ y(tmp);
        yest(tmp)=tmp2*beta;
        plot(x(tmp),yest(tmp),'-','Color',F.loca,'linewidth',3);
    end
    tmp2=[ones(length(x),1) x];
    beta=tmp2 \ y;
    yest=tmp2*beta;
    plot(x,yest,':','Color',F.glob,'linewidth',3);
    set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
    xlim([min(x)-0.2, max(x)]);
    ylim([min(y)-0.2, max(y)]);
else
    set(groot,'defaultAxesColorOrder',map1);
    kmin=2;
    ph=tA(kmin:n,kmin:n)';
    %indPower=find(ph>=(max(max(ph))-0.03));% All scales of 0.03 power diff with max
    % ph(indPower)=2;
    imagesc(ph);
    [k,l]=find(tA==test);
    plot(k-1,l-1,'gs','markerSize',5,'MarkerFaceColor','g')
    
    set(gca,'YDir','normal')
    cmap=map4;
    colormap(cmap)
    % hm=ceil(max(max(ph))*100)/100;
    hm=ceil(prctile(ph(ph<1),99)*100)/100;
    caxis([0 hm])
    h=colorbar('Ticks',[0,hm/2,hm]);%,'location','westoutside');
    set(h,'FontSize',F.fontSize);
    axpos = ax.Position;
    hpos=h.Position;
    hpos(3)=0.5*hpos(3);
    hpos(1)=hpos(1)+0.015;
    h.Position=hpos;
    ax.Position = axpos;
    xlim([1 n-1]);
    ylim([1 n-1]);
    set(gca,'XTick',[2.5,round(n/2)-1,n-1],'YTick',[2.5,round(n/2)-1,n-1],'XTickLabel',[2,round(n/2),n],'YTickLabel',[2,round(n/2),n],'FontSize',F.fontSize);
end
xlabel('X Scales','FontSize',F.fontSize2,...
    'Units', 'normalized','Position', [-0.010, -0.20], 'HorizontalAlignment', 'left');
ylabel('Y Scales','FontSize',F.fontSize2, ...
    'Units', 'normalized', 'Position', [-0.22 -0.02], 'HorizontalAlignment', 'left');

if F.tit
    tit1=strcat('2', F.AB ,'. Multiscale Correlation');
    title([{tit1}; {'Map & Test Statistic'}],'FontSize',F.tfs, ...
        'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left','color','g')
end

set(gca,'FontSize',F.fontSize)
axis('square')
hold off
pos2 = get(ax,'position');
pos2(3:4) = [F.pos(3:4)];
set(ax,'position',pos2);
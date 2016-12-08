function plot_panel3(F,x,y,R2)
ax=subplot('Position',F.pos3);
hold on

n=size(x,1);
test=F.test;
tA=F.tA;
k=F.k;
l=F.l;
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
        plot(x(tmp),yest(tmp),'-','Color',F.loca,'linewidth',2);
    end
    tmp2=[ones(length(x),1) x];
    beta=tmp2 \ y;
    yest=tmp2*beta;
    plot(x,yest,':','Color',F.glob,'linewidth',2);
    set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
    xlim([min(x)-0.2, max(x)]);
    ylim([min(y)-0.2, max(y)+0.1]);
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
if F.type==1
    xlabel('# X Neighbors','FontSize',F.fontSize2+2,...
        'Units', 'normalized','Position', [0.5, -0.08], 'HorizontalAlignment', 'center');
    ylabel('# Y Neighbors','FontSize',F.fontSize2+2, ...
        'Units', 'normalized', 'Position',  [-0.28 0.5], 'HorizontalAlignment', 'center');
end

% if F.tit
%     tit1=strcat('2', F.AB ,'. Multiscale Correlation');
%     title([{tit1}; {'Map & Test Statistic'}],'FontSize',F.tfs, ...
%         'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left','color','g')
%     tit1=strcat('2', F.AB ,'. Multiscale Correlation');
%     tit2=strcat('at (k*, l*) =',{' '},' (',num2str(k),',',{' '},num2str(l), ')');
% end
txt1 = strcat('\color[rgb]{0 1 0} c(MGC)=', num2str(round(100*test)/100));
txt2 = strcat('\color[rgb]{0.5 0.5 0.5} c(Dcorr)=', num2str(round(100*tA(end))/100));
title({txt1;txt2},'FontSize',F.tfs-2,'interpreter','tex');

set(gca,'FontSize',F.fontSize)
axis('square')
hold off
pos2 = get(ax,'position');
pos2(3:4) = [F.pos(3:4)];
set(ax,'position',pos2);
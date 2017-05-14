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
if ~isfield(F,'subplot')
    F.subplot=false;
end
if F.subplot==true
    ax=subplot('Position',F.pos);
else
    ax=figure;
end
hold on
set(groot,'defaultAxesColorOrder',F.map2);
plot(x(:,1),y(:,1),'.','MarkerSize',F.mkSize,'Color',F.gray);
if F.type==1
    if F.subplot==false
        xpos=[0, -0.02];
        ypos=[-0.10 0.5];
    else
        xpos=[0, -0.02];
        ypos=[-0.25 0.5];
    end
    xlabel('Ground Wetness','FontSize',F.fontSize2,...
        'Units', 'normalized','Position', xpos, 'HorizontalAlignment', 'left')
    ylabel('Cloud Density','FontSize',F.fontSize2, ...
        'Units', 'normalized', 'Position',  ypos, 'HorizontalAlignment', 'center')
else
    xlabel('$x$','FontSize',F.fontSize2+5,'Interpreter','latex',...
        'Units', 'normalized','Position', [0.5, -0.02], 'HorizontalAlignment', 'center')
    ylabel('$y$','FontSize',F.fontSize2+5,'Interpreter','latex', ...
        'Units', 'normalized', 'Position',  [-0.02 0.5], 'HorizontalAlignment', 'center')
end

for ind=[1,2,3];
    text(x(F.id(ind),1)+F.hs/3,y(F.id(ind),1)+F.hy(ind)/3,num2str(ind),'fontsize',F.fontSize,'color',F.col)
    plot(x(F.id(ind),1),y(F.id(ind),1),'.','MarkerSize',F.mkSize+3,'Color',F.col);
end

xll=max(x(:,1))-min(x(:,1));
yll=max(y(:,1))-min(y(:,1));
quiver(max(x)+0.25*xll, y(F.id(2),1), 0, y(F.id(3),1)-y(F.id(2),1),'LineWidth',1,'Color',F.col,'MaxHeadSize',0.5,'Autoscale','off');
quiver(max(x)+0.25*xll, y(F.id(3),1), 0, y(F.id(2),1)-y(F.id(3),1),'LineWidth',1,'Color',F.col,'MaxHeadSize',0.5,'Autoscale','off');
quiver(x(F.id(2),1), max(y)+0.25*yll, x(F.id(3),1)-x(F.id(2),1),0,'LineWidth',1,'Color',F.col,'MaxHeadSize',0.5,'Autoscale','off');
quiver(x(F.id(3),1), max(y)+0.25*yll, x(F.id(2),1)-x(F.id(3),1),0,'LineWidth',1,'Color',F.col,'MaxHeadSize',0.5,'Autoscale','off');
txt1 = strcat('$d_{x}(2, 3) =', {' '}, num2str(round(100*abs(x(F.id(3),1)-x(F.id(2),1)))/100),'$');
a=text(x(F.id(2))/2+x(F.id(3))/2,max(y)+0.3*yll,txt1,'VerticalAlignment','bottom','HorizontalAlignment','center','Color',F.col,'Interpreter','latex');
set(a,'FontSize',F.fontSize);
txt1 = strcat('$d_{y}(2, 3) =', {' '}, num2str(round(100*abs(y(F.id(3),1)-y(F.id(2),1)))/100),'$');
a=text(max(x)+0.3*xll,y(F.id(2),1)/2+y(F.id(3),1)/2,txt1,'VerticalAlignment','top','HorizontalAlignment','center','Color',F.col,'Interpreter','latex','rotation',90);
set(a,'FontSize',F.fontSize);

xlim([min(x(:,1))-0.3*xll, max(x(:,1))+0.3*xll]);
ylim([min(y(:,1))-0.3*yll, max(y(:,1))+0.3*yll]);

set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
axis('square')
hold off

% if F.tit
%     tit1=strcat('0', F.AB ,'. Sample Data');
%     title([{tit1}; {' '}], 'Units', 'normalized', ...
%         'Position', [0 1.1], 'HorizontalAlignment', 'left','FontSize',F.tfs);
% end


if F.subplot==false, 
    titletext=CorrSimuTitle(F.type);
    if F.type==1;
        titletext=strcat('A.',{' '}, titletext);
    else
        titletext=strcat('B.',{' '}, titletext);
    end
    h=title(titletext);
    set(h,'FontSize',F.fontSize2+4,'Units', 'normalized' ,'HorizontalAlignment', 'center', 'Units', 'normalized','Position', [0.5 1.1])
    fpath = mfilename('fullpath');
    fpath=strrep(fpath,'\','/');
    findex=strfind(fpath,'/');
    rootDir=fpath(1:findex(end-2));
    pre2=strcat(rootDir,'Figures/');% The folder to save figures
    F.fname=strcat(pre2, 'Fig',num2str(F.type),'Panel',num2str(F.sub));
    print_fig(ax,F);
end
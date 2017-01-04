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
h=subplot('Position',F.pos);
hold all
set(groot,'defaultAxesColorOrder',F.map2);
plot(x,y,'.','MarkerSize',F.mkSize,'Color',F.gray);
if F.type==1
    xlabel('Ground Wetness','FontSize',F.fontSize2,...
        'Units', 'normalized','Position', [0, -0.02], 'HorizontalAlignment', 'left')
    ylabel('Cloud Density','FontSize',F.fontSize2, ...
        'Units', 'normalized', 'Position',  [-0.28 0.5], 'HorizontalAlignment', 'center')
else
    xlabel('$x$','FontSize',F.fontSize2+5,'Interpreter','latex',...
        'Units', 'normalized','Position', [0.5, -0.02], 'HorizontalAlignment', 'center')
    ylabel('$y$','FontSize',F.fontSize2+5,'Interpreter','latex', ...
        'Units', 'normalized', 'Position',  [-0.02 0.5], 'HorizontalAlignment', 'center')
end

for ind=[1,2,3];
    text(x(F.id(ind))+F.hs/3,y(F.id(ind))+F.hy(ind)/3,num2str(ind),'fontsize',F.fontSize,'color',F.col)
    plot(x(F.id(ind)),y(F.id(ind)),'.','MarkerSize',F.mkSize+3,'Color',F.col);
end

xll=max(x)-min(x);
yll=max(y)-min(y);
quiver(max(x)+0.25*xll, y(F.id(2)), 0, y(F.id(3))-y(F.id(2)),'LineWidth',1,'Color',F.col,'MaxHeadSize',0.5,'Autoscale','off');
quiver(max(x)+0.25*xll, y(F.id(3)), 0, y(F.id(2))-y(F.id(3)),'LineWidth',1,'Color',F.col,'MaxHeadSize',0.5,'Autoscale','off');
quiver(x(F.id(2)), max(y)+0.25*yll, x(F.id(3))-x(F.id(2)),0,'LineWidth',1,'Color',F.col,'MaxHeadSize',0.5,'Autoscale','off');
quiver(x(F.id(3)), max(y)+0.25*yll, x(F.id(2))-x(F.id(3)),0,'LineWidth',1,'Color',F.col,'MaxHeadSize',0.5,'Autoscale','off');
txt1 = strcat('$d_{x}(2, 3) =', {' '}, num2str(round(100*abs(x(F.id(3))-x(F.id(2))))/100),'$');
a=text(x(F.id(2))/2+x(F.id(3))/2,max(y)+0.3*yll,txt1,'VerticalAlignment','bottom','HorizontalAlignment','center','Color',F.col,'Interpreter','latex');
set(a,'FontSize',F.fontSize);
txt1 = strcat('$d_{y}(2, 3) =', {' '}, num2str(round(100*abs(y(F.id(3))-y(F.id(2))))/100),'$');
a=text(max(x)+0.3*xll,y(F.id(2))/2+y(F.id(3))/2,txt1,'VerticalAlignment','top','HorizontalAlignment','center','Color',F.col,'Interpreter','latex','rotation',90);
set(a,'FontSize',F.fontSize);

xlim([min(x)-0.3*xll, max(x)+0.3*xll]);
ylim([min(y)-0.3*yll, max(y)+0.3*yll]);

set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
axis('square')

if F.tit
    tit1=strcat('0', F.AB ,'. Sample Data');
    title([{tit1}; {' '}], 'Units', 'normalized', ...
        'Position', [0 1.1], 'HorizontalAlignment', 'left','FontSize',F.tfs);
end


if ~isfield(F,'subprint'), F.subprint=false; end
if F.subprint==true, 
    F.fname=[F.fname, 'a'];
    F.svg=true;
    fpath = mfilename('fullpath');
    fpath=strrep(fpath,'\','/');
    findex=strfind(fpath,'/');
    rootDir=fpath(1:findex(end-2));
    pre2=strcat(rootDir,'Figures/');% The folder to save figures
    F.fname=strcat(pre2, 'Fig',num2str(F.type),'Panel1');
    print_fig(gcf,F)
end
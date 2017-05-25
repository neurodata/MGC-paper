function plot_panel2(F,C,D)

siz=size(C);
n=siz(1);

C2=reshape(C,n^2,1);
D2=reshape(D,n^2,1);

if ~isfield(F,'subplot')
    F.subplot=false;
end
if F.subplot==true
    ax=subplot('Position',F.pos2);
else
    ax=figure;
end
% ax=subplot(s,t,2);
[~,~,RC,RD]=MGCDistTransform(C,D,'mcor');
% RC=(RC<=F.Xmax+1);
% RD=(RD<=F.Ymax+1);

RC=(RC<=F.k);
RD=(RD<=F.l);

ind1=reshape(RC&RD,n^2,1);
hold on
set(groot,'defaultAxesColorOrder',F.map2);



plot(C2,D2,'.','MarkerSize',6,'Color',F.gray);
% if F.sub==3
    plot(C2(ind1==1),D2(ind1==1),'o','MarkerSize',4,'Color',F.loca);
% end



x12=sub2ind([n,n], F.id(1),F.id(2));
x23=sub2ind([n,n], F.id(2),F.id(3));
text(C2(x12)+0.02,D2(x12),'(1, 2)','fontsize',F.fontSize,'color',F.col)
plot(C2(x12),D2(x12),'.','MarkerSize',8,'Color',F.col);

text(C2(x23)+0.02,D2(x23),'(2, 3)','fontsize',F.fontSize,'color',F.col)
plot(C2(x23),D2(x23),'.','MarkerSize',8,'Color',F.col);
hold off
alpha(0.1)

xlim([min(min(C)),1.05]);
ylim([min(min(D)),1.05]);
warning('off','all')
if F.type==1
    if F.subplot==false
        xpos=[0.5, -0.05];
        ypos=[-0.09 0.5];
    else
        xpos=[0.5, -0.05];
        ypos=[-0.21 0.5];
    end
    xlabel('$$d_{x}(x_i,x_j)$$','FontSize',F.fontSize2+2,...
        'Units', 'normalized','Position', xpos, 'HorizontalAlignment', 'center','Interpreter','latex');
    ylabel('$$d_{y}(y_i,y_j)$$','FontSize',F.fontSize2+2, ...
        'Units', 'normalized', 'Position', ypos, 'HorizontalAlignment', 'center','Interpreter','latex');
    set(gca,'XTick',[0,1],'YTick',[0,1],'FontSize',F.fontSize); % Remove x axis tick
else
    set(gca,'XTick',[],'YTick',[],'FontSize',F.fontSize); % Remove x axis tick
end

% if F.sub==2
%     title(strcat('\color[rgb]{0.5 0.5 0.5} c(Dcorr) = ', num2str(round(100*F.tA(end))/100)),'FontSize',F.tfs);
% else
%     txt1 = strcat('(k,l) = (', num2str(F.k),',',num2str(F.l) , ')');
%     txt2 = strcat('c(MGC) = ', num2str(round(100*F.test)/100));
%     title({txt2,txt1},'FontSize',F.tfs,'Color','g'); %,'interpreter','latex');
% %     txt1 = strcat('(k,l) = (', num2str(F.k),',',num2str(F.l) , ')',{', '},'c(MGC) = ', num2str(round(100*F.test)/100));
% %     title(txt1,'FontSize',F.tfs,'Color','g'); %,'interpreter','latex');
% end
txt3 = strcat('\color[rgb]{0.5 0.5 0.5} c(Dcorr) = ', num2str(round(100*F.tA(end))/100));
    title({txt3},'FontSize',F.tfs); %,'interpreter','latex');

set(gca,'FontSize',F.fontSize); % Remove x axis tick
axis('square')
% pos2 = get(ax,'position');
% pos2(3:4) = F.pos(3:4);
% set(ax,'position',pos2);

if F.subplot==false
    fpath = mfilename('fullpath');
    fpath=strrep(fpath,'\','/');
    findex=strfind(fpath,'/');
    rootDir=fpath(1:findex(end-2));
    pre2=strcat(rootDir,'Figures/');% The folder to save figures
    F.fname=strcat(pre2, 'Fig',num2str(F.type),'Panel',num2str(F.sub));
    print_fig(ax,F);
end

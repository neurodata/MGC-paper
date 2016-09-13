% running this code requires running plot_schematic, 
% and breaking after the data and other stuff are loaded

figure(2), clf, 

hold all
set(groot,'defaultAxesColorOrder',map2);
plot(x,y,'.','MarkerSize',mkSize,'Color',gray);
xlabel('x')
ylabel('y')

[I,J]=ind2sub([n,n],find(C_MGC>0.1,1,'first'));
J2=find(mcorrH(J,:)<0,1,'last');
ids=unique([I,J]);
% xx=15;
id=[I,J,J2,J];
id2=[1,2,3,2];
col=[1 .5 0];
hy=[+0.5,-0.5,0];

for ind=[1,2,3]; %length(id)
    hs=0.2;
    text(x(id(ind))+hs,y(id(ind))+hy(ind),num2str(ind),'fontsize',fontSize,'color',col)
    plot(x(id(ind)),y(id(ind)),'.','MarkerSize',mkSize,'Color',col);
end

tname=CorrSimuTitle(type);
findex=strfind(tname,'.');
tname=tname(findex+1:end);
xlim([min(x)-0.2, max(x)]);
ylim([min(y)-0.2, max(y)]);

% title(['0. ', [tname], ' (X,Y)'], 'Units', 'normalized', ...
% title([{'0. Sample Data'}], 'Units', 'normalized', 'Position', [0 1.1], 'HorizontalAlignment', 'left')
set(gca,'XTick',[],'YTick',[],'FontSize',fontSize); % Remove x axis tick
pos=[nan, nan, width, height];
axis('square')



pre2=strcat(rootDir,'docs/images/');% The folder to save figures
F.fname=strcat(pre2, 'spiral');
F.wh=[4 4];
F.PaperPositionMode='auto';
print_fig(gcf,F)

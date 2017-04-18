function []=plot_simulation_visual(opt)
% Author: Cencheng Shen
% CorrVisualPlots()
% CorrVisualPlots(100,2)
% Used to plot figure 0 in the files

%%% File path searching
fpath = mfilename('fullpath');
fpath=strrep(fpath,'\','/');
findex=strfind(fpath,'/');
rootDir=fpath(1:findex(end-2));
addpath(genpath(strcat(rootDir,'Code/')));
pre2=strcat(rootDir,'Figures/Fig');% The folder to save figures

n=30;
dim=1;
if nargin<1
    opt=0;
end

total=20;
figure('units','normalized','position',[0 0 1 1])
s=4;
t=5;
if dim>1
    noise=0;
end

cmap=zeros(3,3);
gr = [0.5,0.5,0.5];
cmap(1,:) = gr;
map1=cmap;
set(groot,'defaultAxesColorOrder',map1);

sz=12;
for type=1:total
    subplot(s,t,type);
    titlechar=strcat(num2str(type),'.',{' '},CorrSimuTitle(type));
    if opt==0
        if type<19
            [x, y]=CorrSampleGenerator(type,n,1,1, noise);
        end
        [x1, y1]=CorrSampleGenerator(type,10*n,dim,1, 0); % Plot 10*n points without noise to highlight the underlying dependency
    end
    if opt==1
        [x, y]=CorrSampleGenerator(type,n,1,1, 1);
    end
    if opt==2
        [x, y]=CorrSampleGenerator(type,n,findDim(type),1, 0);
    end
    sz2=8;
    hold on
    
    if opt==0
        plot(x1(:,1),y1(:,1),'k.','MarkerSize',sz2);
        if type==9
            plot(x1(:,1),y1(:,1),'k.','MarkerSize',sz2*3);
        end
        if type<19
            plot(x(:,1),y(:,1),'.','MarkerSize',sz);
        end
        [a,b]=findRange(type);
        
        xlim(a);
        ylim(b);
        set(gca,'XTick',[]); % Remove x axis ticks
        set(gca,'YTick',[]); % Remove y axis ticks
        set(gca,'box','off','ycolor','w','xcolor','w')
    else
        C=squareform(pdist(x));
%         [~,per]=sort(C(:,1),'ascend');
%         C=C(per,per);
        D=squareform(pdist(y));
%         D=D(per,per);
%         H=eye(n)-ones(n,n)/n;
%         dcov=(H*C*H).*(H*D*H);
%         for i=1:n
%             dcov(i,i)=0;
%         end
%         imagesc(dcov);
        C=C./max(max(C));
        D=D./max(max(D));
        plot(reshape(C,n^2,1),reshape(D,n^2,1),'k.','MarkerSize',5);
        %plot(C(:,1),D(:,1),'k.','MarkerSize',sz2);
        xlim([0,max(max(C))]);
        ylim([0,max(max(D))]);
        if type==1
            set(gca,'XTick',[0,1],'YTick',[0,1]); % Remove x axis ticks
        else
            set(gca,'XTick',[],'YTick',[]); % Remove x axis ticks
        end
    end
    % Specify the axis limit for each type
    
    hold off
    title(titlechar,'FontSize',10)
    axis('square');
    
end
switch opt
    case 0
        h=suptitle('Simulated Example for 20 Dependencies');
    case 1
        h=suptitle('Distance Scatter Plot for 20 Simulations in 1D');
    case 2
        h=suptitle('Distance Scatter Plot for 20 Simulations in HD');
end
set(h,'FontSize',24,'FontWeight','normal');

F.fname=[strcat(pre2, 'SimVisual',num2str(opt))];
F.wh=[8 5]*2;
print_fig(gcf,F)

function [a,b]=findRange(type)
switch type
    case 1
        a=[-1,1];b=[-2.5,2.5];
    case 2
        a=[0,3];b=[-20,40];
    case 3
        a=[-1,1];b=[-250,250];
    case 4
        a=[-3,3];b=[-3,3];
    case 5
        a=[-1,1];b=[-2.5,2.5];
    case 6
        a=[-1,1];b=[-1,2];
    case 7
        a=[-1,1];b=[-1,2];
    case 8
        a=[-6,6];b=[-6,6];
    case 9
        a=[-0.1,1.1];b=[-2,2];
    case 10
        a=[-4,4];b=[-15,10];
    case 11
        a=[-1,1];b=[0,1.5];
    case 12
        a=[-1,1];b=[-2,2];
    case 13
        a=[-1,1];b=[-2,2];
    case 14
        a=[-2,2];b=[-2,2];
    case 15
        a=[-1,1];b=[-1.5,1.5];
    case 16
        a=[-1,1];b=[-1,1];
    case 17
        a=[-6,6];b=[-2,2];
    case 18
        a=[-2,2];b=[-2,2];
    case 19
        a=[-3,3];b=[-3,3];
    case 20
        a=[-3,3];b=[-3,3];
end

function dim=findDim(type)
switch type
    case {1,2,3}
        dim=1000;
    case {5,6,7,15,11,16,17,8}
        dim=20;
    case {4,12,13,19}
        dim=10;
    case {14,18}
        dim=40;
    case {10,9,20}
        dim=100;
end
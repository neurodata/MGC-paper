function []=plot_simulation_permutation(pre1,pre2)

if nargin<1
    pre1='../../Data/Results/'; % The folder to locate data
end
if nargin<2
    pre2='../../Draft/Figures/Fig'; % The folder to save figures
end

%% Set colors
map1=zeros(7,3);
gr = [0,1,0];
ma = [1,0,1];
cy = [0,1,1];
cmap(1,:) = gr;
cmap(2,:) = ma;
cmap(3,:) = cy;
cmap(7,:) = [0 0 0];
dcorr = cmap(1,:);
mcorr = cmap(2,:);
mante = cmap(3,:);
HHG   = [0,0,0];

%
map1(1,:)=dcorr;
map1(2,:)=mcorr; 
map1(3,:)=mcorr; 
map1(4,:)=HHG; 
set(groot,'defaultAxesColorOrder',map1);

figNumber='1DPerm';
filename=strcat(pre1,'CorrSimPermScale1-20Dim1');
load(filename)
figure
x=1:20;
a=2;
if a==1
    ind=[1,2,3,7];
else
    ind=[4,5,6,7];
end
p1=powerP(ind,:);
plot(x,p1(1,:),'.-',x,p1(3,:),'.-',x,p1(2,:),'.:',x,p1(4,:),'.--','LineWidth',2);
legend('Estimated MGC', 'True MGC','Global mcorr','HHG','Location','SouthWest');
set(gca,'FontSize',14);
legend boxoff
xlabel('Function Type','FontSize',15);
ylabel('Testing Power','FontSize',15);
ylim([0,1]);
title('1-Dimensional Simulations at n=60','FontSize',18);
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)

figNumber='HDPerm';
filename=strcat(pre1,'CorrSimPermScale1-20Dim2');
load(filename)
figure
x=1:20;
a=2;
if a==1
    ind=[1,2,3,7];
else
    ind=[4,5,6,7];
end
p1=powerP(ind,:);
plot(x,p1(1,:),'.-',x,p1(3,:),'.-',x,p1(2,:),'.:',x,p1(4,:),'.--','LineWidth',2);
%h=legend('Estimated MGC', 'True MGC','Global Mcorr','Location','SouthWest');
%set(h,'FontSize',12);
set(gca,'FontSize',14);
xlabel('Function Type','FontSize',15);
ylabel('Testing Power','FontSize',15);
ylim([0,1]);
title('High-Dimensional Simulations at n=100','FontSize',18);
F.fname=strcat(pre2, figNumber);
F.wh=[3 2.5]*2;
print_fig(gcf,F)
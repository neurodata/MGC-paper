function [pdsC1, pdsC2, pdsC3]=TibsSimu2Noise(type,n,dim)
%Tibs Simu w.r.t. noise. Compare dCorr, rankdCorr, canonicalCorr
%For y=f(x*A)+eps
%TibsSimu2Noise(1,100, 500);

%parameters
rep=1000; %MC replicates
lim=30; %noise level
display=0; %change to 1 for picture auto-saving
noise=1;

A=ones(dim,1);
%A=A./(ceil(dim*rand(dim,1))); %random decay
for d=1:dim
    A(d)=A(d)/d; %fixed decay
end
dim1=1;

%store intermediate results
pdsC1=zeros(lim+1,1); pdsC2=zeros(lim+1,1);pdsC3=zeros(lim+1,1);%powers for dCorr, rankdCorr, canonical corr
CorA=zeros(rep,lim+1);CorN=zeros(rep,lim+1);
dCor1A=zeros(rep,lim+1);dCor2A=zeros(rep,lim+1);
dCor1N=zeros(rep,lim+1);dCor2N=zeros(rep,lim+1);

%first generate alternative distribution
for l=0:lim
    for r=1:rep
        x=unifrnd(0,1,n,dim);
        xA=x*A;
        
        eps=l/lim*mvnrnd(zeros(n, dim1),eye(dim1)); %Gaussian noise
%         if l~=0
%             eps=trnd(l);
%         else
%             eps=0;
%         end
        
        switch type
            case 1
                y=xA+noise*eps; %linear
            case 2
                y=4*(xA-0.5).^2+10*noise*eps; %Quadratic 1/10
            case 3
                y=128*(xA-1/3).^3+48*(xA-1/3).^2-12*(xA-1/3)+1000*noise*eps; %Cubic 30/1000
            case 4
                y=sin(4*pi*xA/max(xA))+2*noise*eps; %sine 1/2
            case 5
                y=sin(16*pi*xA/max(xA))+noise*eps; %sine 1/8
            case 6
                y=(xA).^0.25+noise*eps; %x^0.25
            case 7
                y=(binornd(1,0.5,n,dim1)*2-1).*(1-(xA/max(xA)*2-1).^2).^0.5+noise/4*eps; %circle
            case 8
                y=(xA>mean(xA))+1*noise*eps; %step 1/5
            case 9
                y=exp(xA)+10*noise*eps; %exponential 1/10
            case 10
                y=log(xA+1)+noise*eps;%log
        end
        x=unifrnd(0,1,n,dim);
        xd=x;
        yd=y;
        CorA(r,l+1)=eucliCorr(xd,yd);
        C=squareform(pdist(xd));
        P=squareform(pdist(yd));
        dCor1A(r,l+1)=distCorr(C,P);
        disA=disToGraph([C P]);
        dCor2A(r,l+1)=distCorr(disA(1:n,1:n), disA(1:n,n+1:2*n));
    end
    
    %then generate null distribution
    for r=1:rep
        x=unifrnd(0,1,n,dim);
        xA=x*A;
        
        eps=l/lim*mvnrnd(zeros(n, dim1),eye(dim1)); %Gaussian noise
%         if l~=0
%             eps=trnd(l);
%         else
%             eps=0;
%         end
        
         switch type
            case 1
                y=xA+noise*eps; %linear
            case 2
                y=4*(xA-0.5).^2+10*noise*eps; %Quadratic 1/10
            case 3
                y=128*(xA-1/3).^3+48*(xA-1/3).^2-12*(xA-1/3)+1000*noise*eps; %Cubic 30/1000
            case 4
                y=sin(4*pi*xA/max(xA))+2*noise*eps; %sine 1/2
            case 5
                y=sin(16*pi*xA/max(xA))+noise*eps; %sine 1/8
            case 6
                y=(xA).^0.25+noise*eps; %x^0.25
            case 7
                y=(binornd(1,0.5,n,dim1)*2-1).*(1-(xA/max(xA)*2-1).^2).^0.5+noise/4*eps; %circle
            case 8
                y=(xA>mean(xA))+1*noise*eps; %step 1/5
            case 9
                y=exp(xA)+10*noise*eps; %exponential 1/10
            case 10
                y=log(xA+1)+noise*eps;%log
        end
        xd=x;
        yd=y;
        CorN(r,l+1)=eucliCorr(xd,yd);
        C=squareform(pdist(xd));
        P=squareform(pdist(yd));
        dCor1N(r,l+1)=distCorr(C,P);
        disA=disToGraph([C P]);
        dCor2N(r,l+1)=distCorr(disA(1:n,1:n), disA(1:n,n+1:2*n));
    end
end

for i=1:lim+1
    dCorT=sort(dCor1A(:,i),'descend');
    cut1=dCorT(ceil(rep*0.05));
    dCorT=sort(dCor2A(:,i),'descend');
    cut2=dCorT(ceil(rep*0.05));
    dCorT=sort(CorA(:,i),'descend');
    cut3=dCorT(ceil(rep*0.05));
    pdsC1(i)=mean(dCor1N(:,i)>cut1);
    pdsC2(i)=mean(dCor2N(:,i)>cut2);
    pdsC3(i)=mean(CorN(:,i)>cut3);
end

if size(type,1)>1
    type=100;
end
filename=strcat('TibsSimuInd2NoiseType',num2str(type),'N',num2str(n));
save(filename,'pdsC1','pdsC2','pdsC3','type','n','noise','rep','lim','dim','dCor1N','dCor2N','CorN','dCor1A','dCor2A','CorA');

%display and save picture. By default do not display.
if display~=0
    filename=strcat('TibsSimuInd2NoiseType',num2str(type),'N',num2str(n));
    titlechar='Data';
    switch type
        case 1
            titlechar='Linear';
        case 2
            titlechar='Quadratic';
        case 3
            titlechar='Cubic';
        case 4
            titlechar='Sine Period 1/2';
        case 5
            titlechar='Sine Period 1/8';
        case 6
            titlechar='X\^(1/4)';
        case 7
            titlechar='Circle';
        case 8
            titlechar='Step Function';
        case 9
            titlechar='Exponential';
        case 10
            titlechar='Log';
    end
    
    %figure
    xaxis=1:lim+1;
    %%%Plot with canonical correlation. But pdsC3 has zero power for hd data when d > n
    %plot(xaxis,pdsC1,'b*-',xaxis,pdsC2,'ro-',xaxis,pdsC3,'g.-','LineWidth',1);
    %legend('Distance Correlation','Ranking Distance Correlation','Canonical Correlation', 'Location','NorthEast');
    plot(xaxis, pdsC1,'b*-',xaxis,pdsC2,'ro-','LineWidth',1);
    legend('Distance Correlation','Ranking Distance Correlation','Location','NorthEast');
    
    %figure title/labels
    ylabel('Testing Power \beta');
    xlim([1 lim+1]);
    ylim([0 1]);
    xlabel('Noise Level');
    ax=gca;
    ax.XTick=1:lim/10:lim+1;
    ax.XTickLabel=0:noise/10:noise;
    %xlabel('Degree of Freedom');
    titleStr = strcat(titlechar,' Type Independence Test with Noise at n=', num2str(n), ' and dim= ',num2str(dim));
    title(titleStr);
    saveas(gcf,filename,'jpeg');
end

function disA=disToGraph(dis) %transform from distance to ranking
n=size(dis,1);
[~,ind]=sort(dis(1:n,1:n));
[~,ind2]=sort(dis(1:n,n+1:2*n));
disA=zeros(n,2*n);
for i=1:n
    v=ind(2:n,i);
    w=ind2(2:n,i);
    for j=1:n-1
        disA(v(j),i)=j;
        %disA(i,v(j),K)=max(disA(i,v(j),K),disA(v(j),i,K));
        disA(w(j),n+i)=j;
        %disA(i,n+w(j),K)=max(disA(i,n+w(j),K),disA(w(j),n+i,K));
    end
end

function corr = distCorr(X,Y) %calculate dCorr
X=(X+X')/2;
Y=(Y+Y')/2;
n=size(X,1);
H=eye(n)-(1/n)*ones(n,n);
X=H*X*H;
Y=H*Y*H;
corr=abs(sum(sum(X.*(Y)))); %need absolute for ranking DCorr, but not necessary for usual DCorr
corr=sqrt(corr/(norm(X,'fro')*norm(Y,'fro')));

function corr = eucliCorr(X,Y) %not really useful...yet...
[~,~,corr]=canoncorr(X,Y);
corr=corr(1);

% function corr = graphCorr(X,Y)
% corr=sum(sum(X.*(Y)));
% corr=sqrt(corr/(norm(X,'fro')*norm(Y,'fro')));
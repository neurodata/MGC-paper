rm(list = ls())
library(ggplot2)
library(reshape2)
library(lattice)
library(R.matlab)
source("CorrSimTitles.R")
total=20;
optionA=2;

#load("CorrPermDistTestTypeData.RData")

#Fig1
for (j in (1:total)){
  filename = paste("CorrIndTestType", j,"N100Dim1.mat",sep = "");
  a=readMat(filename);
  plotName=paste("pl", j,sep = "");
  
  titleName=CorrSimTitles(a$type);n=a$n;numRange=a$numRange;p4=a$power4;
  p1L=rep(0,a$lim);p2L=rep(0,a$lim);p3L=rep(0,a$lim);
  p1M=rep(0,a$lim);p2M=rep(0,a$lim);p3M=rep(0,a$lim);
  for (i in (1:a$lim)){
    p1M[i]=max(a$power1[2:numRange[i],2:numRange[i],i]);
    p2M[i]=max(a$power2[2:numRange[i],2:numRange[i],i]);
    p3M[i]=max(a$power3[2:numRange[i],2:numRange[i],i]);
    p1L[i]=max(a$power1[numRange[i],numRange[i],i]);
    p2L[i]=max(a$power2[numRange[i],numRange[i],i]);
    p3L[i]=max(a$power3[numRange[i],numRange[i],i]);
  }
  if (optionA==1){
    plotData=data.frame(LGCmcorr=p1M,mcorr=p1L,dcorr=p2L,Mantel=p3L,HHG=t(p4),numRange=t(numRange));
  } else {
    plotData=data.frame(LGCmcorr=p1M,LGCdcorr=p2M,LGCMantel=p3M,mcorr=p1L,dcorr=p2L,Mantel=p3L,HHG=t(p4),numRange=t(numRange));
  }
  plotData <- melt(plotData,id="numRange");
  colnames(plotData) <- c("numRange","Method","Power");
  
  if (optionA==1){
    lineText=c(rep("solid", 1), rep("dotted", 4));
    colorText=c("red","red","blue","cyan","green");
    legendText=c("LGC by mcorr","mcorr","dcorr","Mantel","HHG");
  } else {
    lineText=c(rep("solid", 3), rep("dotted", 4));
    colorText=c("red","blue","cyan","red","blue","cyan","green");
    legendText=c("LGC by mcorr","LGC by dcorr","LGC by Mantel","mcorr","dcorr","Mantel","HHG")
  }
  pl=ggplot(data=plotData,aes(x=numRange,y=Power,color=Method,linetype=Method))+geom_line(size=2)+geom_point(size=1)+ylab("P-Value")+ggtitle(titleName)+scale_linetype_manual(values=lineText)+scale_colour_manual(values=colorText)+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),plot.title = element_text(size=16,face="bold"),legend.title = element_text(size=16, face="bold"),legend.text = element_text(size=12),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ ylim(0, 1);
  assign(plotName,pl);
}
source("multiplot.R")
#par(oma = c(4, 1, 1, 1))
multiplot(pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9, pl10, pl11, pl12, pl13, pl14, pl15, pl16, pl17, pl18, pl19, pl20, cols=4);
#par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#legend("bottom", inset = 0,legend=legendText, col=colorText, lty=lineText, cex=1, horiz = TRUE)

#Fig3
thres=0.8;
interval=0.01;total=20;
xaxis=seq(0,1,interval);
profile=matrix(0,7,length(xaxis));
for (j in (1:total)){
  filename = paste("CorrIndTestType", j,"N100Dim1.mat",sep = "");
  a=readMat(filename);
  
  n=a$n;numRange=a$numRange;p4=a$power4;
  p1L=rep(0,a$lim);p2L=rep(0,a$lim);p3L=rep(0,a$lim);
  p1M=rep(0,a$lim);p2M=rep(0,a$lim);p3M=rep(0,a$lim);
  for (i in (1:a$lim)){
    p1M[i]=max(a$power1[2:numRange[i],2:numRange[i],i]);
    p2M[i]=max(a$power2[2:numRange[i],2:numRange[i],i]);
    p3M[i]=max(a$power3[2:numRange[i],2:numRange[i],i]);
    p1L[i]=max(a$power1[numRange[i],numRange[i],i]);
    p2L[i]=max(a$power2[numRange[i],numRange[i],i]);
    p3L[i]=max(a$power3[numRange[i],numRange[i],i]);
  }
  ind=c(which(p1M>=thres),which(p2L>=thres),which(p3L>=thres),which(p4>=thres),a$lim);
  pos=min(ind);
  if (optionA==1){
    power=c(p1M[pos],p1L[pos],p2L[pos],p3L[pos],p4[pos]);
  } else {
    power=c(p1M[pos],p1L[pos],p2L[pos],p3L[pos],p4[pos],p2M[pos],p3M[pos]);
  }
  pmax=max(power);
  for (k in (1:length(power))){
    tmp=(pmax-power[k]);
    tmpInd=ceiling(tmp/interval)+1;
    profile[k,tmpInd:(1/interval+1)]=profile[k,tmpInd:(1/interval+1)]+1;
  }
}
profile=profile/total;
sumP=ceiling(rowMeans(profile)*1000)/1000;
if (optionA==1){
  plotData=data.frame(LGCmcorr=profile[1,],mcorr=profile[2,],dcorr=profile[3,],Mantel=profile[4,],HHG=profile[5,],numRange=xaxis);
} else {
  plotData=data.frame(LGCmcorr=profile[1,],LGCdcorr=profile[6,],LGCMantel=profile[7,],mcorr=profile[2,],dcorr=profile[3,],Mantel=profile[4,],HHG=profile[5,],numRange=xaxis);
}
plotData <- melt(plotData,id="numRange");
colnames(plotData) <- c("Interval","Method","Profiles");
if (optionA==1){
  lineText=c(rep("solid", 1), rep("dotted", 4));
  colorText=c("red","red","blue","cyan","green");
  legendText=c("LGC by mcorr","mcorr","dcorr","Mantel","HHG");
} else {
  lineText=c(rep("solid", 3), rep("dotted", 4));
  colorText=c("red","blue","cyan","red","blue","cyan","green");
  legendText=c("LGC by mcorr","LGC by dcorr","LGC by Mantel","mcorr","dcorr","Mantel","HHG")
}
titleName="Performance Profiles for Dimension 1";
pl=ggplot(data=plotData,aes(x=Interval,y=Profiles,color=Method,linetype=Method))+geom_line(size=2)+geom_point(size=1)+xlab("Difference with the Best Method")+ylab("Relative Performance")+ggtitle(titleName)+scale_linetype_manual(values=lineText)+scale_colour_manual(values=colorText)+theme(legend.position="right",plot.title = element_text(size=16,face="bold"),legend.title = element_text(size=16, face="bold"),legend.text = element_text(size=12),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ ylim(0, 1);
pl
#800*600

#Fig4
interval=0.01;total=20;
xaxis=seq(0,1,interval);
limN=10;sumP=matrix(0,7,limN);
for (ll in (1:limN)){
  profile=matrix(0,7,length(xaxis));
  thres=ll/limN;
  for (j in (1:total)){
    filename = paste("CorrIndTestType", j,"N100Dim1.mat",sep = "");
    a=readMat(filename);
    n=a$n;numRange=a$numRange;p4=a$power4;
    p1L=rep(0,a$lim);p2L=rep(0,a$lim);p3L=rep(0,a$lim);
    p1M=rep(0,a$lim);p2M=rep(0,a$lim);p3M=rep(0,a$lim);
    for (i in (1:a$lim)){
      p1M[i]=max(a$power1[2:numRange[i],2:numRange[i],i]);
      p2M[i]=max(a$power2[2:numRange[i],2:numRange[i],i]);
      p3M[i]=max(a$power3[2:numRange[i],2:numRange[i],i]);
      p1L[i]=max(a$power1[numRange[i],numRange[i],i]);
      p2L[i]=max(a$power2[numRange[i],numRange[i],i]);
      p3L[i]=max(a$power3[numRange[i],numRange[i],i]);
    }
    ind=c(which(p1M>=thres),which(p2L>=thres),which(p3L>=thres),which(p4>=thres),a$lim);
    pos=min(ind);
    if (optionA==1){
      power=c(p1M[pos],p1L[pos],p2L[pos],p3L[pos],p4[pos]);
    } else {
      power=c(p1M[pos],p1L[pos],p2L[pos],p3L[pos],p4[pos],p2M[pos],p3M[pos]);
    }
    pmax=max(power);
    for (k in (1:length(power))){
      tmp=(pmax-power[k]);
      tmpInd=ceiling(tmp/interval)+1;
      profile[k,tmpInd:(1/interval+1)]=profile[k,tmpInd:(1/interval+1)]+1;
    }
  }
  profile=profile/total;
  sumP[,ll]=ceiling(rowMeans(profile)*1000)/1000;
}
xaxis=seq(1/limN,1,1/limN);
if (optionA==1){
  plotData=data.frame(LGCmcorr=sumP[1,],mcorr=sumP[2,],dcorr=sumP[3,],Mantel=sumP[4,],HHG=sumP[5,],numRange=xaxis);
} else {
  plotData=data.frame(LGCmcorr=sumP[1,],LGCdcorr=sumP[6,],LGCMantel=sumP[7,],mcorr=sumP[2,],dcorr=sumP[3,],Mantel=sumP[4,],HHG=sumP[5,],numRange=xaxis);
}
plotData <- melt(plotData,id="numRange");
colnames(plotData) <- c("Interval","Method","Profiles");
if (optionA==1){
  lineText=c(rep("solid", 1), rep("dotted", 4));
  colorText=c("red","red","blue","cyan","green");
  legendText=c("LGC by mcorr","mcorr","dcorr","Mantel","HHG");
} else {
  lineText=c(rep("solid", 3), rep("dotted", 4));
  colorText=c("red","blue","cyan","red","blue","cyan","green");
  legendText=c("LGC by mcorr","LGC by dcorr","LGC by Mantel","mcorr","dcorr","Mantel","HHG")
}
titleName="Area Under Curve of Performance Profiles for Dimension 1";
pl=ggplot(data=plotData,aes(x=Interval,y=Profiles,color=Method,linetype=Method))+geom_line(size=2)+geom_point(size=1)+xlab("Threshold of Power")+ylab("Area Under Curve")+ggtitle(titleName)+scale_linetype_manual(values=lineText)+scale_colour_manual(values=colorText)+theme(legend.position="right",plot.title = element_text(size=16,face="bold"),legend.title = element_text(size=16, face="bold"),legend.text = element_text(size=12),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ ylim(0, 1);
pl



#Fig5
rm(list = ls())
total=20;
optionA=1;
library(ggplot2)
library(reshape2)
source("CorrSimTitles.R")
library(R.matlab)
for (j in (1:total)){
  filename = paste("CorrIndTestDimType", j,"N100Dim.mat",sep = "");
  a=readMat(filename);
  plotName=paste("pl", j,sep = "");
  
  titleName=CorrSimTitles(a$type);n=a$n;numRange=a$dimRange;
  p1L=rep(0,a$lim);p2L=rep(0,a$lim);p3L=rep(0,a$lim);
  p1M=rep(0,a$lim);p2M=rep(0,a$lim);p3M=rep(0,a$lim);p4=a$power4;
  for (i in (1:a$lim)){
    p1M[i]=max(a$power1[2:n,2:n,i]);
    p2M[i]=max(a$power2[2:n,2:n,i]);
    p3M[i]=max(a$power3[2:n,2:n,i]);
    p1L[i]=max(a$power1[n,n,i]);
    p2L[i]=max(a$power2[n,n,i]);
    p3L[i]=max(a$power3[n,n,i]);
  }
  if (optionA==1){
    plotData=data.frame(LGCmcorr=p1M,mcorr=p1L,dcorr=p2L,Mantel=p3L,HHG=t(p4),numRange=t(numRange));
  } else {
    plotData=data.frame(LGCmcorr=p1M,LGCdcorr=p2M,LGCMantel=p3M,mcorr=p1L,dcorr=p2L,Mantel=p3L,HHG=t(p4),numRange=t(numRange));
  }
  plotData <- melt(plotData,id="numRange");
  colnames(plotData) <- c("numRange","Method","Power");
  
  if (optionA==1){
    lineText=c(rep("solid", 1), rep("dotted", 4));
    colorText=c("red","red","blue","cyan","green");
    legendText=c("LGC by mcorr","mcorr","dcorr","Mantel","HHG");
  } else {
    lineText=c(rep("solid", 3), rep("dotted", 4));
    colorText=c("red","blue","cyan","red","blue","cyan","green");
    legendText=c("LGC by mcorr","LGC by dcorr","LGC by Mantel","mcorr","dcorr","Mantel","HHG")
  }
  pl=ggplot(data=plotData,aes(x=numRange,y=Power,color=Method,linetype=Method))+geom_line(size=2)+geom_point(size=1)+ylab("P-Value")+ggtitle(titleName)+scale_linetype_manual(values=lineText)+scale_colour_manual(values=colorText)+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank(),plot.title = element_text(size=16,face="bold"),legend.title = element_text(size=16, face="bold"),legend.text = element_text(size=12),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ ylim(0, 1);
  assign(plotName,pl);
}
source("multiplot.R")
#par(oma = c(4, 1, 1, 1))
multiplot(pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9, pl10, pl11, pl12, pl13, pl14, pl15, pl16, pl17, pl18, pl19, pl20, cols=4);
#par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#legend("bottom", inset = 0,legend=legendText, col=colorText, lty=lineText, cex=1, horiz = TRUE)


#Fig7
thres=0.5;
interval=0.01;total=20;
xaxis=seq(0,1,interval);
profile=matrix(0,7,length(xaxis));
for (j in (1:total)){
  filename = paste("CorrIndTestDimType", j,"N100Dim.mat",sep = "");
  a=readMat(filename);
  
  n=a$n;numRange=a$dimRange;p4=a$power4;
  p1L=rep(0,a$lim);p2L=rep(0,a$lim);p3L=rep(0,a$lim);
  p1M=rep(0,a$lim);p2M=rep(0,a$lim);p3M=rep(0,a$lim);
  for (i in (1:a$lim)){
    p1M[i]=max(a$power1[2:n,2:n,i]);
    p2M[i]=max(a$power2[2:n,2:n,i]);
    p3M[i]=max(a$power3[2:n,2:n,i]);
    p1L[i]=max(a$power1[n,n,i]);
    p2L[i]=max(a$power2[n,n,i]);
    p3L[i]=max(a$power3[n,n,i]);
  }
  ind=c(which(p1M>=thres),which(p2L>=thres),which(p3L>=thres),which(p4>=thres),1);
  pos=max(ind);
  if (optionA==1){
    power=c(p1M[pos],p1L[pos],p2L[pos],p3L[pos],p4[pos]);
  } else {
    power=c(p1M[pos],p1L[pos],p2L[pos],p3L[pos],p4[pos],p2M[pos],p3M[pos]);
  }
  pmax=max(power);
  for (k in (1:length(power))){
    tmp=(pmax-power[k]);
    tmpInd=ceiling(tmp/interval)+1;
    profile[k,tmpInd:(1/interval+1)]=profile[k,tmpInd:(1/interval+1)]+1;
  }
}
profile=profile/total;
sumP=ceiling(rowMeans(profile)*1000)/1000;
if (optionA==1){
  plotData=data.frame(LGCmcorr=profile[1,],mcorr=profile[2,],dcorr=profile[3,],Mantel=profile[4,],HHG=profile[5,],numRange=xaxis);
} else {
  plotData=data.frame(LGCmcorr=profile[1,],LGCdcorr=profile[6,],LGCMantel=profile[7,],mcorr=profile[2,],dcorr=profile[3,],Mantel=profile[4,],HHG=profile[5,],numRange=xaxis);
}
plotData <- melt(plotData,id="numRange");
colnames(plotData) <- c("Interval","Method","Profiles");
if (optionA==1){
  lineText=c(rep("solid", 1), rep("dotted", 4));
  colorText=c("red","red","blue","cyan","green");
  legendText=c("LGC by mcorr","mcorr","dcorr","Mantel","HHG");
} else {
  lineText=c(rep("solid", 3), rep("dotted", 4));
  colorText=c("red","blue","cyan","red","blue","cyan","green");
  legendText=c("LGC by mcorr","LGC by dcorr","LGC by Mantel","mcorr","dcorr","Mantel","HHG")
}
titleName="Performance Profiles for Increasing Dimension";
pl=ggplot(data=plotData,aes(x=Interval,y=Profiles,color=Method,linetype=Method))+geom_line(size=2)+geom_point(size=1)+xlab("Difference with the Best Method")+ylab("Relative Performance")+ggtitle(titleName)+scale_linetype_manual(values=lineText)+scale_colour_manual(values=colorText)+theme(legend.position="right",plot.title = element_text(size=16,face="bold"),legend.title = element_text(size=16, face="bold"),legend.text = element_text(size=12),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ ylim(0, 1);
pl

#Fig8
interval=0.01;total=20;
xaxis=seq(0,1,interval);
limN=10;sumP=matrix(0,7,limN);
for (ll in (1:limN)){
  profile=matrix(0,7,length(xaxis));
  thres=ll/limN;
  for (j in (1:total)){
    filename = paste("CorrIndTestDimType", j,"N100Dim.mat",sep = "");
    a=readMat(filename);
    n=a$n;numRange=a$dimRange;p4=a$power4;
    p1L=rep(0,a$lim);p2L=rep(0,a$lim);p3L=rep(0,a$lim);
    p1M=rep(0,a$lim);p2M=rep(0,a$lim);p3M=rep(0,a$lim);
    for (i in (1:a$lim)){
      p1M[i]=max(a$power1[2:n,2:n,i]);
      p2M[i]=max(a$power2[2:n,2:n,i]);
      p3M[i]=max(a$power3[2:n,2:n,i]);
      p1L[i]=max(a$power1[n,n,i]);
      p2L[i]=max(a$power2[n,n,i]);
      p3L[i]=max(a$power3[n,n,i]);
    }
    ind=c(which(p1M>=thres),which(p2L>=thres),which(p3L>=thres),which(p4>=thres),1);
    pos=max(ind);
    if (optionA==1){
      power=c(p1M[pos],p1L[pos],p2L[pos],p3L[pos],p4[pos]);
    } else {
      power=c(p1M[pos],p1L[pos],p2L[pos],p3L[pos],p4[pos],p2M[pos],p3M[pos]);
    }
    pmax=max(power);
    for (k in (1:length(power))){
      tmp=(pmax-power[k]);
      tmpInd=ceiling(tmp/interval)+1;
      profile[k,tmpInd:(1/interval+1)]=profile[k,tmpInd:(1/interval+1)]+1;
    }
  }
  profile=profile/total;
  sumP[,ll]=ceiling(rowMeans(profile)*1000)/1000;
}
xaxis=seq(1/limN,1,1/limN);
if (optionA==1){
  plotData=data.frame(LGCmcorr=sumP[1,],mcorr=sumP[2,],dcorr=sumP[3,],Mantel=sumP[4,],HHG=sumP[5,],numRange=xaxis);
} else {
  plotData=data.frame(LGCmcorr=sumP[1,],LGCdcorr=sumP[6,],LGCMantel=sumP[7,],mcorr=sumP[2,],dcorr=sumP[3,],Mantel=sumP[4,],HHG=sumP[5,],numRange=xaxis);
}
plotData <- melt(plotData,id="numRange");
colnames(plotData) <- c("Interval","Method","Profiles");
if (optionA==1){
  lineText=c(rep("solid", 1), rep("dotted", 4));
  colorText=c("red","red","blue","cyan","green");
  legendText=c("LGC by mcorr","mcorr","dcorr","Mantel","HHG");
} else {
  lineText=c(rep("solid", 3), rep("dotted", 4));
  colorText=c("red","blue","cyan","red","blue","cyan","green");
  legendText=c("LGC by mcorr","LGC by dcorr","LGC by Mantel","mcorr","dcorr","Mantel","HHG")
}
titleName="Area Under Curve of Performance Profiles for Increasing Dimension";
pl=ggplot(data=plotData,aes(x=Interval,y=Profiles,color=Method,linetype=Method))+geom_line(size=2)+geom_point(size=1)+xlab("Threshold of Power")+ylab("Area Under Curve")+ggtitle(titleName)+scale_linetype_manual(values=lineText)+scale_colour_manual(values=colorText)+theme(legend.position="right",plot.title = element_text(size=16,face="bold"),legend.title = element_text(size=16, face="bold"),legend.text = element_text(size=12),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ ylim(0, 1);
pl




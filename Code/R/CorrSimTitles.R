CorrSimTitles <- function(type){ # Calculate local dCorr
  
  if (type==1){
    titlechar=' 1.Linear';
  }
  if (type==2){
    titlechar=' 2.Cubic';
  }
  if (type==3){
    titlechar=' 3.Step Function';
  }
  if (type==4){
    titlechar=' 4.Exponential';
  }
  if (type==5){
    titlechar=' 5.Joint Normal';
  }
  if (type==6){
    titlechar=' 6.Quadratic';
  }
  if (type==7){
    titlechar=' 7.W Shape';
  }
  if (type==8){
    titlechar=' 8.Two Parabolas';
  }
  if (type==9){
    titlechar=' 9.Fourth Root';
  }
  if (type==10){
    titlechar=' 10.Logarithmic';
  }
  if (type==11){
    titlechar=' 11.Circle';
  }
  if (type==12){
    titlechar=' 12.Ellipse';
  }
  if (type==13){
    titlechar=' 13.Spiral';
  }
  if (type==14){
    titlechar=' 14.Square';
  }
  if (type==15){
    titlechar=' 15.Diamond';
  }
  if (type==16){
    titlechar=' 16.Sine Period 1/2';
  }
  if (type==17){
    titlechar=' 17.Sine Period 1/8';
  }
  if (type==18){
    titlechar=' 18.Multiplicative Noise';
  }
  if (type==19){
    titlechar=' 19.Uncorrelated Binomial';
  }
  if (type==20){
    titlechar=' 20.Independent Clouds';
  }
  
  return(titlechar);
}
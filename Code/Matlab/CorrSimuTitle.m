function titlechar=CorrSimuTitle(type)
%%%Used to call the character names of each dependency in the simulation
switch type
    case 1
        titlechar=' 1.Linear';
    case 2
        titlechar=' 2.Cubic';
    case 3
        titlechar=' 3.Step Function';
    case 4
        titlechar=' 4.Exponential';
    case 5
        titlechar=' 5.Joint Normal';
    case 6
        titlechar=' 6.Quadratic';
    case 7
        titlechar=' 7.W Shape';
    case 8
        titlechar=' 8.Two Parabolas';
    case 9
        titlechar=' 9.Fourth Root';
    case 10
        titlechar=' 10.Logarithmic';
    case 11
        titlechar=' 11.Circle';
    case 12
        titlechar=' 12.Ellipse';
    case 13
        titlechar=' 13.Spiral';
    case 14
        titlechar=' 14.Square';
    case 15
        titlechar=' 15.Diamond';
    case 16
        titlechar=' 16.Sine Period 1/2';
    case 17
        titlechar=' 17.Sine Period 1/8';
    case 18
        titlechar=' 18.Multiplicative Noise';
    case 19
        titlechar=' 19.Uncorrelated Binomial';
    case 20
        titlechar=' 20.Independent Clouds';
end
function titlechar=CorrSimuTitle(type)
%%%Used to call the character names of each dependency in the simulation
switch type
    case 1
        titlechar='Linear';
    case 2
        titlechar='Exponential';
    case 3
        titlechar='Cubic';
    case 4
        titlechar='Joint Normal';
    case 5
        titlechar='Step Function';
    case 6
        titlechar='Quadratic';
    case 7
        titlechar='W Shape';
    case 8
        titlechar='Spiral';
    case 9
        titlechar='Bernoulli';
    case 10
        titlechar='Logarithmic';
    case 11
        titlechar='Fourth Root';
    case 12
        titlechar='Sine Period 4\pi';
    case 13
        titlechar='Sine Period 16\pi';
    case 14
        titlechar='Square';
    case 15
        titlechar='Two Parabolas';
    case 16
        titlechar='Circle';
    case 17
        titlechar='Ellipse';
    case 18
        titlechar='Diamond';
    case 19
        titlechar='Multiplicative';
    case 20
        titlechar='Independent';
end
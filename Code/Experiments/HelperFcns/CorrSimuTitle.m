function titlechar=CorrSimuTitle(type)
%%%Used to call the character names of each dependency in the simulation
switch type
    case 1
        titlechar='1. Linear';
    case 2
        titlechar='2. Exponential';
    case 3
        titlechar='3. Cubic';
    case 4
        titlechar='4. Joint Normal';
    case 5
        titlechar='5. Step Function';
    case 6
        titlechar='6. Quadratic';
    case 7
        titlechar='7. W Shape';
    case 8
        titlechar='8. Spiral';
    case 9
        titlechar='9. Bernoulli';
    case 10
        titlechar='10. Logarithmic';
    case 11
        titlechar='11. Fourth Root';
    case 12
        titlechar='12. Sine Period 4\pi';
    case 13
        titlechar='13. Sine Period 16\pi';
    case 14
        titlechar='14. Square';
    case 15
        titlechar='15. Two Parabolas';
    case 16
        titlechar='16. Circle';
    case 17
        titlechar='17. Ellipse';
    case 18
        titlechar='18. Diamond';
    case 19
        titlechar='19. Multiplicative';
    case 20
        titlechar='20. Independent';
end
function [legPlastic, bodyPlastic] = definePlasticPropsTypHex()
%Properties for the plastic of the leg segments, taken from the solidworks
%models
legPlastic.length{1} = [24; 41.3; 100; 100]/1000; %m
legPlastic.length{2} = legPlastic.length{1};
legPlastic.length{3} = legPlastic.length{1};
legPlastic.length{4} = legPlastic.length{1};
legPlastic.length{5} = [17.4+24; 21.2; 80; 80]/1000;
legPlastic.length{6} = legPlastic.length{5};

legPlastic.width{1} = [50; 25.4; 25.4; 25.4]/1000;
legPlastic.width{2} = legPlastic.width{1};
legPlastic.width{3} = legPlastic.width{1};
legPlastic.width{4} = legPlastic.width{1};
legPlastic.width{5} = [50; 25.4; 25.4; 25.4]/1000;
legPlastic.width{6} = legPlastic.width{5};

legPlastic.height{1} = [23.8; 50.8; 25.4; 25.4]/1000;
legPlastic.height{2} = legPlastic.height{1};
legPlastic.height{3} = legPlastic.height{1};
legPlastic.height{4} = legPlastic.height{1};
legPlastic.height{5} = [23.8; 50.8; 25.4; 25.4]/1000;
legPlastic.height{6} = legPlastic.height{5};

legPlastic.COM{1} = [0.81, 2.60, 4.87, 3.15;...
    0.00, 1.26, 0.00, 0.00;...
    0.00, 0.00, 0.00, 0.00]/100;
legPlastic.COM{2} = legPlastic.COM{1};
legPlastic.COM{3} = legPlastic.COM{1};
legPlastic.COM{4} = legPlastic.COM{1};
legPlastic.COM{5} = [1.22, 0, 0;...
    1.05, 1.27, 0;...
    3.9, 0, 0;...
    2.3, 0, 0]'/100;
legPlastic.COM{6} = legPlastic.COM{5};

legPlastic.mass{1} = [6.4; 18.7; 29.7; 25.2]/1000;
legPlastic.mass{2} = legPlastic.mass{1};
legPlastic.mass{3} = legPlastic.mass{1};
legPlastic.mass{4} = legPlastic.mass{1};
legPlastic.mass{5} = [3.9+6.4; 12.7; 23.6; 20.7]/1000; %kg
legPlastic.mass{6} = legPlastic.mass{5};

%Properties for the plastic of the thorax, taken from the solidworks models
bodyPlastic.mass = [37.86; 92.19; 68.46]/1000;
bodyPlastic.length = [50; 121.9; 107]/1000; %m
bodyPlastic.width = [170; 170; 175]/1000;
bodyPlastic.height = [6; 75; 6]/1000;
bodyPlastic.COM = [26.2, 0, 0;...
    64.5, 39.8, 0;...
    48.7, 0, 0]'/1000; %COM considered from the most posterior face
end
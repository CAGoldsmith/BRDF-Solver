%Clarissa Goldsmith
%West Virginia University

%Assuming we're using the MX-28T as the actuators
actuator.mass = 82.7/1000; %kg - MX-28T
actuator.length = 50.6/1000; %m
actuator.width = 35.5/1000; %m
actuator.height = 35.6/1000; %m
actuator.hornWidth = 40/1000; %m
actuator.hornCOMLength = 16.1/1000; %m
actuator.I = 0.00943;

%Segment masses are:
plasM = [4.6; 3.9; 6.4; 12.7; 23.6; 20.7]/1000; %kg
segM = [actuator.mass+plasM(1); actuator.mass+plasM(2); plasM(3); actuator.mass*2 + plasM(4); plasM(5); plasM(6)]; %kg
function I_leg = findIMatrixGeneric(robot,comPoints)
numLegs = robot.numLegs;
for i=1:numLegs %For each leg...
    numBodies = robot.legObj{i}.numBodies-2;
    segM = robot.legObj{i}.mass(2:end-1);
    segW = robot.legObj{i}.width(2:end-1);
    segH = robot.legObj{i}.height(2:end-1);
    segL = robot.legObj{i}.lengths(2:end-1);
    for b=1:numBodies
        I_leg{i}{b} = zeros(3,3);
        I_leg{i}{b}(1,1) = 1/12 * segM(b) * (segW(b)^2 + segH(b)^2); %Ixx
        I_leg{i}{b}(2,2) = 1/12 * segM(b) * (segL(b)^2 + segW(b)^2); %Iyy
        I_leg{i}{b}(3,3) = 1/12 * segM(b) * (segL(b)^2 + segH(b)^2); %Izz
    end  
end
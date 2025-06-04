function [I_leg,I_body] = findIMatrix(robot)
numLegs = robot.numLegs; 

for i=1:numLegs %For each leg...


    numBodies = robot.legObj{i}.numBodies;
    segL = robot.legObj{i}.lengths;
    plasL = robot.legObj{i}.plasticL;
    plasW = robot.legObj{i}.plasticW;
    plasH = robot.legObj{i}.plasticH;
    plasM = robot.legObj{i}.plasticMass;
    plasCOM = robot.legObj{i}.plasticCOM;
    segM = robot.legObj{i}.mass;
    actuator = robot.actuator;

    %% Calculate the CoM of each segment with reference to the most proximal joint it attaches to

    if i == 5 || i == 6
        j = 2;
        k = 1;
        %ThC Bridge:
        x1 = actuator.hornWidth/2; y1 = actuator.hornCOMLength;
        x2 = plasCOM(1,1);
        x3 = x2 + actuator.length/2;
        x1_s = -(y1*cosd(38)+(x1+y1*tand(38)))*actuator.mass;
        y1_s = (x1+y1*tand(38))*cosd(38)*actuator.mass;
        x2_s = x2*sind(38)*plasM(2);
        y2_s = x2*cosd(38)*plasM(2);
        x3_s = x3*sind(38)*actuator.mass;
        y3_s = x3*cosd(38)*actuator.mass;
        COMjr(:,1) = [(x2_s + x3_s)/segM(1);...
            (y2_s - y3_s)/segM(1);...
            0];

        %Coxa:
        COMjr(:,2) = [(plasCOM(1,2)*plasM(2) + (plasCOM(1,2)+actuator.length/2)*(135/1000))/segM(2);...
            0; 0];

    else
        j = 1;
        k = 0;
        %Coxa:
        COMjr(:,1) = [segL(1+j)-actuator.length/2; 0; 0];
    end

    %Trochanter:
    COMjr(:,2+k) = plasCOM(:,2+k);

    %Femur:
    COMjr(1,3+k) = (actuator.hornWidth/2*actuator.mass + (actuator.hornWidth + plasCOM(1,3+k))*plasM(3+k) +...
        (actuator.hornWidth + plasL(3+k) + actuator.height/2)*actuator.mass)/segM(3+j);
    COMjr(2,3+k) = (-actuator.length/2*actuator.mass + plasCOM(2,3+k)*plasM(3+k) + actuator.length/2*actuator.mass)/segM(3+j);
    COMjr(3,3+k) = 0;

    %Tibia
    COMjr(:,4+k) = plasCOM(:,4+k);

    %Tarsus
    COMjr(:,5+k) = plasCOM(:,5+k);

    %% Compute the moments/products of inertia for each segment
    if i == 5 || i == 6
        %THC BRIDGE:
        %First, need to get the non-rotated moments and products
        I_bridgeNR(1,1) = 1/12 * plasM(1) * (plasW(1)^2 + plasH(1)^2) +... %I about CoM for plastic
            1/12 * actuator.mass * (actuator.hornWidth^2 + actuator.height^2) +... %I about CoM for lateral actuator
            actuator.I; %Contribution of actuator's internal I

        I_bridgeNR(2,2) = 1/12 * plasM(1) * (plasL(1)^2 + plasW(1)^2) +... %I about CoM for plastic
            plasM(1) * abs(plasCOM(1,1))^2 +... %parallel axis thm for plastic
            1/12 * actuator.mass * (actuator.length^2 + actuator.height^2) +... %I about CoM for lateral actuator
            actuator.mass * abs(plasCOM(1,1)+actuator.length/2)^2 -... %parallel axis thm for lateral actuator
            segM(1) * abs((plasM(1)*plasCOM(1,1)+actuator.mass*(plasCOM(1,1)+actuator.length/2))/segM(1))^2; %parallel axis thm for whole segment

        I_bridgeNR(3,3) = 1/12 * plasM(1) * (plasL(1)^2 + plasH(1)^2) +... %I about CoM for plastic
            plasM(1) * abs(plasCOM(1,1))^2 +... %parallel axis thm for plastic
            1/12 * actuator.mass * (actuator.length^2 + actuator.hornWidth^2) +... %I about CoM for lateral actuator
            actuator.mass * abs(plasCOM(1,1)+actuator.length/2)^2 -... %parallel axis thm for lateral actuator
            segM(1) * abs((plasM(1)*plasCOM(1,1)+actuator.mass*(plasCOM(1,1)+actuator.length/2))/segM(1))^2; %parallel axis thm for whole segment

        %Because the only deviation from the CoM is in the x direction, the non-rotated products will all be zero

        %Now, do the rotation, using code taken from oneLegInvKin.m
        w = [0;0;1];
        wHat = hat(w);
        v = [0;0;0];
        exp_w = @(theta,wHat) eye(3)+wHat*sind(theta)+wHat^2*(1-cosd(theta));
        g = @(theta,wHat,v) [exp_w(theta,wHat), (eye(3)-exp_w(theta,wHat))*wHat*v; 0 0 0 1];
        g_bridge = g(52,wHat,v);
        g_bridge = g_bridge(1:3,1:3);

        I_bridge = g_bridge*I_bridgeNR*g_bridge';

        %COXA:
        %Moments of inertia
        I_coxa(1,1) = 1/12 * plasM(2) * (plasW(2)^2 + plasH(2)^2) +... %I around CoM for plastic
            1/12 * (135/1000) * (actuator.height^2 + actuator.hornWidth^2)+... %I around CoM for actuator
            actuator.I; %Moment contribution from servo
        %No parallel axis thm here because CoM in line with joint point for y and z
        %coords
        I_coxa(2,2) = 1/12 * plasM(2) * (plasW(2)^2 + plasL(2)^2) +... %I around CoM for plastic
            plasM(2) * abs(plasCOM(1,2))^2 +... %Parallel axis thm for plastic
            1/12 * (135/1000) * (actuator.length^2 + actuator.hornWidth^2) +... %I around CoM or actuator
            (135/1000) * abs(plasCOM(1,2) + actuator.length/2)^2 -... %parallel axis thm for plastic
            segM(2) * abs(COMjr(1,2))^2; %Parallel axis thm for whole segment

        I_coxa(3,3) = 1/12 * plasM(2) * (plasH(2)^2 + plasL(2)^2) +... %I around CoM for plastic
            plasM(2) * abs(plasCOM(1,2))^2 +... %Parallel axis thm for plastic
            1/12 * (135/1000) * (actuator.height^2 + actuator.length^2) +... %I around CoM for actuator
            (135/1000) * abs(plasCOM(1,2) + actuator.length/2)^2 -... %Parallel axis thm for actuator
            segM(2) * abs(COMjr(1,2))^2; %Parallel axis thm for whole segment

        %Then, the products of inertia
        %All products of inertia for the coxa are zero, because the only deviation
        %is in x
    else
        I_bridge = 0;

        %COXA:
        %Start with the moments of inertia
        I_coxa(1,1) = 1/12 * (135/1000) * (actuator.height^2 + actuator.hornWidth^2); %Ixx
        I_coxa(2,2) = 1/12 * (135/1000) * (actuator.length^2 + actuator.hornWidth^2); %Iyy
        I_coxa(3,3) = 1/12 * (135/1000) * (actuator.height^2 + actuator.length^2); %Izz
        %Then, the products of inertia
        %All products of inertia for the coxa are zero
    end

    %TROCHANTER:
    %Moments of inertia:
    I_troc(1,1) = 1/12 * plasM(2+k) * (plasW(2+k)^2 + plasH(2+k)^2);
    I_troc(2,2) = 1/12 * plasM(2+k) * (plasW(2+k)^2 + plasL(2+k)^2);
    I_troc(3,3) = 1/12 * plasM(2+k) * (plasL(2+k)^2 + plasH(2+k)^2) + actuator.I;
    %Product of inertia:
    %All products of inertia for the trochanter are zero

    %FEMUR
    %This one is the bear
    %Moments of inertia:
    I_fem(1,1) = (1/12 * actuator.mass * (actuator.height^2 + actuator.length^2) +... %I about CoM for medial actuator
        actuator.mass * norm([-actuator.length/2; 0])^2 +... %parallel axis thm for medial actuator
        1/12 * plasM(3+k) * (plasW(3+k)^2 + plasH(3+k)^2) +... %I abt CoM for plastic
        plasM(3+k) * norm([plasCOM(2,3+k);0])^2 +... %parallel axis thm for plastic
        1/12 * actuator.mass * (actuator.hornWidth^2 + actuator.length^2) +...%I abt CoM for lateral actuator
        actuator.mass * norm([actuator.length/2; 0])^2) -...%parallel axis thm for lateral actuator
        (segM(3+j) * norm(COMjr(2:3,3+k))^2)+...%parallel axis thm for whole segment
        actuator.I; %TrF rotates around the x axis, so add in servo I

    I_fem(2,2) = (1/12 * actuator.mass * (actuator.hornWidth^2 + actuator.height^2) +... %I about CoM for medial actuator
        actuator.mass * norm([actuator.hornWidth/2; 0])^2 +... %parallel axis thm for medial actuator
        1/12 * plasM(3+k) * (plasW(3+k)^2 + plasL(3+k)^2)^2 +... %I abt CoM for plastic
        plasM(3+k) * norm([actuator.hornWidth + plasCOM(1,3+k); 0])^2 +... %parallel axis thm for plastic
        1/12 * actuator.mass * (actuator.hornWidth^2 + actuator.height^2) +... %I abt CoM for lateral actuator
        actuator.mass * norm([actuator.hornWidth + plasL(3+k) + actuator.height/2; 0])^2) -... %parallel axis thm for lateral actuator
        (segM(3+j)* norm([COMjr(1,3+k); COMjr(3,3+k)])^2); %parallel axis thm for whole segment

    I_fem(3,3) = (1/12 * actuator.mass * (actuator.length^2 + actuator.hornWidth^2) +... %I about CoM for medial actuator
        actuator.mass * norm([actuator.hornWidth/2; -actuator.length/2])^2 +... %parallel axis thm for medial actuator
        1/12 * plasM(3+k) * (plasL(3+k)^2 + plasH(3+k)^2) +... %I abt CoM for plastic
        plasM(3+k) * norm([actuator.hornWidth + plasCOM(1,3+k); plasCOM(2,3+k)])^2 +... %parallel axis thm for plastic
        1/12 * actuator.mass * (actuator.length^2 + actuator.height^2) +... %I abt CoM for lateral actuator
        actuator.mass * norm([actuator.hornWidth + plasL(3+k) + actuator.height/2; actuator.length/2])^2) -... %parallel axis thm for lateral actuator
        segM(3+j) * norm(COMjr(1:2,3+k))^2; %parallel axis thm for whole segment

    %Products of inertia:
    I_fem(2,1) = -(actuator.mass * -(segL(3+j)/2 - actuator.hornWidth/2) * -actuator.length/2 +...
        segM(3+j) * (actuator.hornWidth + 2.04/100 - segL(3+j)/2) * plasCOM(2,3+k) +...
        actuator.mass * (segL(3+j)/2 - actuator.hornWidth/2) * actuator.length/2);
    I_fem(1,2) = I_fem(2,1);

    %Other products of inertia should be zero, as the z coord of the CoM
    %doesn't leave the x-axis for any of the portions of the segment

    %TIBIA
    %Moments of Inertia:
    I_tib(1,1) = 1/12 * segM(4+j) * (plasW(4+k)^2 + plasH(4+k)^2);
    I_tib(2,2) = 1/12 * segM(4+j) * (plasL(4+k)^2 + plasW(4+k)^2);
    I_tib(3,3) = 1/12 * segM(4+j) * (plasL(4+k)^2 + plasH(4+k)^2) + actuator.I;

    %Products of inertia: All are zero

    %TARSUS
    %Moments of inertia:
    I_tar(1,1) = 1/12 * segM(5+j) * (plasW(5+k)^2 + plasH(5+k)^2);
    I_tar(2,2) = 1/12 * segM(5+j) * (plasL(5+k)^2 + plasW(5+k)^2);
    I_tar(3,3) = 1/12 * segM(5+j) * (plasL(5+k)^2 + plasH(5+k)^2);
    %Products of inertia: All are zero

    if i == 5 || i == 6
        I_leg{i} = {I_bridge; I_coxa; I_troc; I_fem; I_tib; I_tar};
    else
        I_leg{i} = {I_coxa; I_troc; I_fem; I_tib; I_tar};
    end
end

%% Find the center of mass for the thorax
bodySegMass = robot.plastic.mass;
bodySegL = robot.plastic.length;
bodySegW = robot.plastic.width;
bodySegH = robot.plastic.height;
bodySegLocCOM = robot.plastic.COM; 
bodySegGlobCOM = bodySegLocCOM + [sum(bodySegL(2:3)),bodySegH(2),0; bodySegL(3),0,0; 0,0,0]';

bodyCOM(1,1) = (bodySegMass(1)*bodySegGlobCOM(1,1) + bodySegMass(2)*bodySegGlobCOM(1,2) + bodySegMass(3)*bodySegGlobCOM(1,3))/sum(bodySegMass);
bodyCOM(2,1) = (bodySegMass(1)*bodySegGlobCOM(1,2) + bodySegMass(2)*bodySegGlobCOM(2,2)) / sum(bodySegMass);
bodyCOM(3,1) = 0;

%% Find the vectors between the thorax COM and the COM of every thorax segment
rG1G = bodySegGlobCOM(:,1)-bodyCOM;
rG2G = bodySegGlobCOM(:,2)-bodyCOM;
rG3G = bodySegGlobCOM(:,3)-bodyCOM;

%% Find the thorax moments/products of inertia
%Moments of Inertia
I_body(1,1) = 1/12*bodySegMass(1)*(bodySegH(1)^2+bodySegW(1)^2) +... %Moment of segment 1 around its COM
    + bodySegMass(1)*norm(rG1G(2:3))^2 +... %Parallel axis thm for segment 1 around body COM
    + 1/12*bodySegMass(2)*(bodySegH(1)^2+bodySegW(2)^2) +... %Moment of segment 2 around its COM
    + bodySegMass(2)*norm(rG2G(2:3))^2 +... %Parallel axis thm for segment 2 around body COM
    + 1/12*bodySegMass(3)*(bodySegH(3)^2+bodySegW(3)^2) +... %Moment of segment 3 around its COM
    + bodySegMass(3)*norm(rG3G(2:3))^2; %Parallel axis thm for segment 3 around body COM

I_body(2,2) = 1/12*bodySegMass(1)*(bodySegL(1)^2+bodySegW(1)^2) +... %Moment of segment 1 around its COM
    + bodySegMass(1)*norm(rG1G(1))^2 +... %Parallel axis thm for segment 1 around body COM
    + 1/12*bodySegMass(2)*(bodySegL(2)^2+bodySegW(2)^2) +... %Moment of segment 2 around its COM
    + bodySegMass(2)*norm(rG2G(1))^2 +... %Parallel axis thm for segment 2 around body COM
    + 1/12*bodySegMass(3)*(bodySegL(3)^2+bodySegW(3)^2) +... %Moment of segment 3 around its COM
    + bodySegMass(3)*norm(rG3G(1))^2; %Parallel axis thm for segment 3 around body COM

I_body(3,3) = 1/12*bodySegMass(1)*(bodySegH(1)^2+bodySegL(1)^2) +... %Moment of segment 1 around its COM
    + bodySegMass(1)*norm(rG1G(1:2))^2 +... %Parallel axis thm for segment 1 around body COM
    + 1/12*bodySegMass(2)*(bodySegH(1)^2+bodySegL(2)^2) +... %Moment of segment 2 around its COM
    + bodySegMass(2)*norm(rG2G(1:2))^2 +... %Parallel axis thm for segment 2 around body COM
    + 1/12*bodySegMass(3)*(bodySegH(3)^2+bodySegL(3)^2) +... %Moment of segment 3 around its COM
    + bodySegMass(3)*norm(rG3G(1:2))^2; %Parallel axis thm for segment 3 around body COM

%Products of inertia
I_body(1,2) = -(bodySegMass(1)*rG1G(1)*rG1G(2) +  bodySegMass(2)*rG2G(1)*rG2G(2) + bodySegMass(3)*rG3G(1)*rG3G(2));
I_body(2,1) = I_body(1,2);
%The rest of the products will be zero

end

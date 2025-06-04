function [] = design_drosophibot(savePathBase,configName,saveName,swingDuration,stanceDuty,stepHeight,floorLevel,springK,maxJointVelocity,bodyTranslation,bodyRotation,IMatrixFunName,actuatorName,terrainShape,contralateralPhase,ipsalateralPhase,runNeuralParams)
load([savePathBase '\variables\configs\' configName ' Config.mat'])
%%
mirrored = [];
% swingDuration = 0.5; %s
maxStanceDuration = 6; %s
% stanceDuty = 0.6; %Percentage of time stance takes in a step
% stepHeight = [0.07;0.07;0.06;.06;.07;.07];
% stepProfile = 'simp';
stepProfile = 'biom';

% if physicalRobot
%     floorLevel = [-.175;-.175;-.19;-.19;-.165;-.165];
% %         floorLevel = [-.175;-.175;-.19;-.19;-.175;-.175];
% else
%     floorLevel = -.175*ones(6,1);
% end

% bodyTranslation = 8.5e-2; %m
% bodyRotation = 0.2; %rad
% maxJointVelocity = 1.5; %3; %rad/s
constantAngleJoints = [];
isRobot = true;%,run
loopParameters = [];
tune = tic;
% springK = 0.307; %Nm/rad
% springK = 0;

%Positive x: forward
%Positive y: upward
%Positive z: right

%Leg 1: LH
%Leg 2: RH
%Leg 3: LM
%Leg 4: RM
%Leg 5: LF
%Leg 6: RF

swingPeriodPlot = ['T_swing = ' num2str(swingDuration)];
springKPlot = ['k_spring = ' num2str(springK)];
swingPeriodPlotTitle = ['T_{swing} = ' num2str(swingDuration)];
springKPlotTitle = ['k_{spring} = ' num2str(springK)];

if ~exist([savePathBase '\' saveName '\' terrainShape '\' springKPlot '\' swingPeriodPlot '\phi_I = ' num2str(ipsalateralPhase) '\phi_C = ' num2str(contralateralPhase)])
    mkdir([savePathBase '\' saveName '\' terrainShape '\' springKPlot '\' swingPeriodPlot '\phi_I = ' num2str(ipsalateralPhase) '\phi_C = ' num2str(contralateralPhase)])
end
savePath = [savePathBase '\' saveName '\' terrainShape '\' springKPlot '\' swingPeriodPlot '\phi_I = ' num2str(ipsalateralPhase) '\phi_C = ' num2str(contralateralPhase)];

%
load([savePathBase '\variables\actuators\' actuatorName '.mat'])

% %Properties for the MX-28
% MX28.mass = 82.7/1000; %kg - MX-28T
% % MX28.mass = 135/1000; %kg - MX-64T
% MX28.length = 50.6/1000; %m
% MX28.width = 35.5/1000; %m
% MX28.height = 35.6/1000; %m
% MX28.hornWidth = 40/1000; %m
% MX28.hornCOMLength = 16.1/1000; %m
% MX28.c = 4.5;
% MX28.k = 8.4;
% MX28.I = 0.00943;

% Pull the plastic properties from the appropriate function
if contains(configName, '2')
    %For Drosophibot 2:
    [legPlastic, bodyPlastic] = definePlasticPropsDroso2();
elseif contains(configName, 'Hex')
    %For typical hexapod configuration
    [legPlastic, bodyPlastic] = definePlasticPropsTypHex();
else
    legPlastic = [];
    bodyPlastic = [];
end

%This function scrapes all of the necessary info to recreate the organism
%structure from the Animatlab file into a kinematicOrganism variable:
robot = kinematicOrganism(proj_file,organism_name,bodies,joints,ctr,...
    swingDuration,maxStanceDuration,bodyTranslation,bodyRotation,constantAngleJoints,...
    maxJointVelocity,mirrored,isRobot,true,loopParameters,legPlastic,bodyPlastic,actuator);

%Compute the Jacobians for each leg segment
[~,~,~] = robot.computeAllJacobians([],true);
%close all


%%
boundingFactor = 0;
postureScaleFactor = robot.legObj{robot.numLegs/2}.lengths(end-1);
ballScaleFactor = mean([abs(robot.pJoints{1}(1,1)-robot.pJoints{5}(1,1)), robot.legObj{robot.numLegs/2}.lengths(end-1)]); 
%If terrain is a ball, define the ball radius and the position of the
%center in the world frame
if contains(terrainShape,'Ball')
    ballParams.r = 0.003/6.5375e-04*ballScaleFactor;
    if contains(configName, 'Drosophibot')
        ballParams.center = [robot.pJoints{3}(1,1);-0.0037./6.5375e-04*ballScaleFactor+.07; 0];
    elseif contains(configName, 'Scale')
        ballParams.center = [robot.pJoints{3}(1,1)/2;-0.0037./6.5375e-04*ballScaleFactor; 0];
    end
else
    ballParams = [];
end
[anglesRest,footPointsRest] = robot.findRestingPosture(jointVecs,floorLevel,boundingFactor,postureScaleFactor,terrainShape,ballParams);
robot.setLowLevelParams();
% configMat = cell2mat(configs);


signMat = [];
for i=1:numLegs
    if rem(i,2) == 0
        signMat(1,i) = 1;
        signMat(5,i) = 1;
    else
        signMat(1,i) = -1;
        signMat(5,i) = -1;
    end
end
signMat(2,:) = -ones(1,numLegs);
signMat(3,:) = ones(1,numLegs);
signMat(4,:) = zeros(1,numLegs);

saveas(gcf,[savePathBase '\' configName '\' terrainShape '\Resting Posture.fig']);
% input('Press enter to close all figures and continue.\n')
close all
%%
% pBase = [-.3;floorLevel(1);0];
pBase = [footPointsRest{1}(1,1);floorLevel(1);0];
[rotationMat,translationMat,kinCell,stKin,~,~,~,pepConfigs,aepConfigs] = robot.walkingInverseKinematics(jointVecs,pBase,stanceDuty,stepHeight,stepProfile,terrainShape,ballParams);
close all
if ~exist('jointLimits')
    [br,bc] = size(bodies);
    jointLimits = cell(br,bc);
end

if contains(runNeuralParams,'Yes')
    neuralParamsSavePath = [savePathBase '\' saveName '\' terrainShape];
    robot.computeAllRobotAdapterGains(jointLimits,neuralParamsSavePath);
    close all
    robot.setMotionAngles(jointVecs,pepConfigs,aepConfigs,rotationMat,translationMat,neuralParamsSavePath);
end

for i=[1 4 6] %CHANGE AT SOME POINT IF DOING MORE THAN SIX LEGS
    figure
    for j=1:robot.legObj{i}.numBodies-2
        if contains(terrainShape,'Flat')
            plot(kinCell{i,7,8}(j,:),'-o')
        else
            plot(kinCell{i,2}(j,:),'-o')
        end
        hold on
        xlim([0 length(kinCell{1,1,1})])
    end
    title(['Straight Walking Inverse Kinematics' ' Leg ' num2str(i)])
end

%% NOW THAT WE HAVE THE IK, WE CAN COMPUTE THE JOINT TORQUES, INTERNAL FORCES, AND GROUND REACTION FORCES

numTranslations = length(unique(translationMat));
numRotations = length(unique(rotationMat));

sampsPerStep = size(kinCell{1,1,1},2);
stepPeriod = (swingDuration/(1-stanceDuty));
dt = stepPeriod/sampsPerStep;

%One option: just look at specific walking directions (this is straight
%forward)
if contains(terrainShape,'Flat')
    transVec = 7;
    rotVec = 8;
else
    transVec = 2;
    rotVec = 1;
end
%Another option: look at multiple/all walking directions (this is all of
%them)
% transVec = 1:numTranslations;
% rotVec = 1:numRotations;

stepSize = 10;
toPlot = 1;
legNames = {'LH','RH','LM','RM','LF','RF'};

%Calculate the total mass of the organism's legs based on the masses from the
%animatlab file
orgLegMassTot = 0;
for i=1:numLegs
    orgLegMassTot = orgLegMassTot + sum(robot.legObj{i}.mass(2:end-1));
end

for tran=transVec
    savePathDir = [savePath '\Trans = ' num2str(tran)];
    if ~exist(savePathDir)
        mkdir(savePathDir)
    end
    for rot=rotVec %FOR A GIVEN WALKING DIRECTION:

        for i = 1:numLegs
            stanceIDsPt1 = 1:1:length(stKin{i,tran,rot})/2;
            stanceIDsPt2 = fliplr(length(kinCell{2,tran,rot}):-1:length(kinCell{2,tran,rot})-(length(stKin{i,tran,rot})/2-1));
            stanceIDs{i} = [stanceIDsPt1 stanceIDsPt2];
        end

        %Shift the legs an appropriate amount to create the desired gait,
        %as determined by ispalateralPhase and contralateralPhase

        %Leg 1 (LH) will always be considered at phase 0
        for L=1:numLegs
            if L == 1
                phaseAmnt(L) = 0;
            else
                if ~rem(L,2)
                    phaseAmnt(L) = phaseAmnt(L-1) + contralateralPhase;
                else
                    phaseAmnt(L) = phaseAmnt(L-2) + ipsalateralPhase;
                end
            end
            if phaseAmnt(L) > 1
                phaseAmnt(L) = phaseAmnt(L) - 1;
            end
            shiftAmnt(L) = round(phaseAmnt(L)*sampsPerStep);
            if phaseAmnt(L) > 0
                kinCell{L,tran,rot} = circshift(kinCell{L,tran,rot},shiftAmnt(L),2);
                stanceIDs{L} = stanceIDs{L} + shiftAmnt(L);
                for qq = 1:length(stanceIDs{L})
                    if stanceIDs{L}(qq) > sampsPerStep
                        stanceIDs{L}(qq) = stanceIDs{L}(qq) - sampsPerStep;
                    end
                end
            end
        end

        timeVector = 1:1:sampsPerStep;

        F1(round(sampsPerStep/stepSize)) = struct('cdata',[],'colormap',[]);
        F2(round(sampsPerStep/stepSize)) = struct('cdata',[],'colormap',[]);        

        nn = 0;
        for n=timeVector
            %If you want to plot the legs during each frame, specify a
            %figure handle:
            if toPlot
                if rem(n,stepSize) == 1
                    figHandle1 = figure;
                    title(sprintf('Frame %i',n))
                    nn = nn + 1;
                else
                    figHandle1 = [];
                end
            else
                figHandle1 = [];
            end

            for i=1:numLegs %FOR EACH LEG:
                %Extract the theta values into a new array
                theta{n}{i} = kinCell{i,tran,rot}(:,n);
                %Thetas stored from most proximal to least proximal
                if toPlot
                    if rem(n,stepSize) == 1
                        robot.legObj{i}.setToPlot(true);
                    else
                        robot.legObj{i}.setToPlot(false);
                    end
                else
                    robot.legObj{i}.setToPlot(false);
                end

                %Calculate the foot, joint, and COM points for this frame
                [jacobianBodyFrame{n}{i},footPoint{n}{i},jointPoints{n}{i},comPoints{n}{i},omegas{n}{i},~,g_b{n}] = robot.legObj{i}.computeJacobian(theta{n}{i},figHandle1);
                
                for j=1:length(robot.legObj{i}.joints)-2
                    jointPointsPlot{i}{j}(:,n) = jointPoints{n}{i}(:,j);
                    footPointPlot{i}(:,n) = footPoint{n}{i};
                end
                thetaPlot{i}(n,:) = theta{n}{i}';
            end

            %Calculate where the center of mass should be based on the CoM
            %of all the points and the masses
            CoMLocNum = zeros(3,1);
            for i=1:numLegs
                    c = comPoints{n}{i}(:,2:end-1);
                    massRange = 2:robot.legObj{i}.numBodies-1;
                for j = massRange
                    c(:,j-1) = c(:,j-1)*robot.legObj{i}.mass(j);
                end
            CoMLocNum = CoMLocNum + sum(c,2);
            end
            CoMCoord(:,n) = CoMLocNum/(orgLegMassTot+robot.legObj{1}.mass(1));

            %For each leg, call the findMassChars function(s) for the legs to calculate the
            %inertia matrices for each leg
            %Find the I matrices for all of the legs and the thorax
            if n == 1
                if strcmp(IMatrixFunName,'findIMatrix')
                    %For Drosophibot 2:
                    [I_leg,~] = findIMatrix(robot);
                elseif strcmp(IMatrixFunName,'findIMatrixGiga')
                    %For Giga Drosophila:
                    I_leg = findIMatrixGiga(robot,comPoints);
                elseif strcmp(IMatrixFunName,'findIMatrixTypHex')
                    %For a typical hexapod:
                    I_leg = findIMatrixTypHex(robot);
                else
                    I_leg = findIMatrixGeneric(robot,comPoints);
                end
            end

            if rem(n,stepSize) == 1 %Create plots of the side profile and overhead profile at this timestep, to be stitched together in a gif
                ylim([footPointsRest{3}(3)*1.25 footPointsRest{4}(3)*1.25])
                zlim([min(floorLevel)*1.1 abs(min(floorLevel))/2])
                xlim([footPointsRest{1}(1,1)*1.5 footPointsRest{numLegs}(1,1)*1.5])
                hold on
                plot3(CoMCoord(1,n),-CoMCoord(3,n),CoMCoord(2,n),'diamond')
                if contains(terrainShape,'Ball')
                    [sphereX,sphereY,sphereZ] = sphere;
                    sphereX = sphereX*ballParams.r;
                    sphereY = sphereY*ballParams.r;
                    sphereZ = sphereZ*ballParams.r;
                    surf(sphereX+ballParams.center(1),-(sphereZ+ballParams.center(3)),sphereY+ballParams.center(2));
                end
                F1(nn) = getframe(figHandle1);
                view(0,90)
                plot3([footPoint{n}{2}(1) footPoint{n}{3}(1),footPoint{n}{6}(1) footPoint{n}{2}(1)],[-footPoint{n}{2}(3) -footPoint{n}{3}(3),-footPoint{n}{6}(3) -footPoint{n}{2}(3)],[footPoint{n}{2}(2) footPoint{n}{3}(2),footPoint{n}{6}(2) footPoint{n}{2}(2)],'-o')
                F2(nn) = getframe(figHandle1);
                if ~exist([savePathDir '\Walking Frames'])
                    mkdir([savePathDir '\Walking Frames'])
                end
                saveas(figHandle1,[savePathDir '\Walking Frames\Frame ' num2str(n)])
                close(figHandle1)
            end
        end

        %Convert the figures into a gif to check the walking pattern
        if toPlot
            for f=1:sampsPerStep/stepSize
                im1{f} = frame2im(F1(f));
                im2{f} = frame2im(F2(f));
                [Aim1,map1] = rgb2ind(im1{f},256);
                [Aim2,map2] = rgb2ind(im2{f},256);
                if f == 1
                    imwrite(Aim1,map1,[savePathDir '\robotWalkingPatternSide.gif'],'gif','LoopCount',Inf,'DelayTime',1/10);
                    imwrite(Aim2,map2,[savePathDir '\robotWalkingPatternTop.gif'],'gif','LoopCount',Inf,'DelayTime',1/10);
                else
                    imwrite(Aim1,map1,[savePathDir '\robotWalkingPatternSide.gif'],'gif','WriteMode','append','DelayTime',1/10);
                    imwrite(Aim2,map2,[savePathDir '\robotWalkingPatternTop.gif'],'gif','WriteMode','append','DelayTime',1/10);
                end
            end
        end
        %Calculate the the omegas for each joint and velocities of each
        %body by averaging over the previous and future timesteps
        for n=timeVector
            for i=1:numLegs
                segNum = robot.legObj{i}.numBodies-2;
                if n > 1 && n < sampsPerStep
                    thetaDscalar{i}(:,n) = (kinCell{i,tran,rot}(:,n+1)-kinCell{i,tran,rot}(:,n-1))/(2*dt);
                    segLinVel{n}{i} = (comPoints{n+1}{i}(:,2:end-1) - comPoints{n-1}{i}(:,2:end-1))/(2*dt);
                elseif n == 1
                    thetaDscalar{i}(:,n) = (kinCell{i,tran,rot}(:,n+1)-kinCell{i,tran,rot}(:,end-(n-1)))/(2*dt);
                    segLinVel{n}{i} = (comPoints{n+1}{i}(:,2:end-1) - comPoints{end-(n-1)}{i}(:,2:end-1))/(2*dt);
                elseif n == sampsPerStep
                    thetaDscalar{i}(:,n) = (kinCell{i,tran,rot}(:,1)-kinCell{i,tran,rot}(:,n-1))/(2*dt);
                    segLinVel{n}{i} = (comPoints{1}{i}(:,2:end-1) - comPoints{n-1}{i}(:,2:end-1))/(2*dt);
                end
                for j=1:segNum
                    thetaD{n}{i}(:,j) = thetaDscalar{i}(j,n) * omegas{n}{i}(:,j);
                end
                for j=1:segNum
                    segAngVel{n}{i}(:,j) = sum(thetaD{n}{i}(:,1:j),2);
                end
            end
        end

        %Now repeat cycling through the timeVector, this time calculating
        %acceleration
        for n=timeVector
            for i=1:numLegs
                segNum = robot.legObj{i}.numBodies-2;
                if contains(configName, 'Massless')
                    segLinAcc{n}{i} = zeros(3,segNum);
                    segAngAcc{n}{i} = zeros(3,segNum);
                else
                    if n > 1 && n < sampsPerStep
                        thetaDDscalar{i}(:,n) = (thetaDscalar{i}(:,n+1)-thetaDscalar{i}(:,n-1))/(2*dt);
                        segLinAcc{n}{i} = (segLinVel{n+1}{i}-segLinVel{n-1}{i})/(2*dt);
                    elseif n == 1
                        thetaDDscalar{i}(:,n) = (thetaDscalar{i}(:,n+1)-thetaDscalar{i}(:,end))/(2*dt);
                        segLinAcc{n}{i} = (segLinVel{n+1}{i}-segLinVel{end}{i})/(2*dt);
                    else
                        thetaDDscalar{i}(:,n) = (thetaDscalar{i}(:,1)-thetaDscalar{i}(:,n-1))/(2*dt);
                        segLinAcc{n}{i} = (segLinVel{1}{i}-segLinVel{n-1}{i})/(2*dt);
                    end
                    for j=1:segNum
                        segAngAcc{n}{i}(:,j) = thetaDDscalar{i}(j,n) * omegas{n}{i}(:,j);
                        for jj = 1:j
                            segAngAcc{n}{i}(:,j) = segAngAcc{n}{i}(:,j) + cross(thetaD{n}{i}(:,jj),omegas{n}{i}(:,j));
                        end
                    end
                end
            end
        end

        %Find the minimum theta for each joint to serve as the eq. pos for
        %the springs
        for i=1:numLegs
            segNum = robot.legObj{i}.numBodies-2;  
            for j=1:segNum
                if contains(robot.legObj{i}.joints{j+1},'ThC3') 
                    if contains(robot.legObj{i}.joints{j+1},'LF') ||  contains(robot.legObj{i}.joints{j+1},'RF')
                        eqAngle{i}(j,1) = anglesRest{i}(j);
                    elseif contains(robot.legObj{i}.joints{j+1},'LM') ||  contains(robot.legObj{i}.joints{j+1},'RM')
                        eqAngle{i}(j,1) = 0;
                    else
                        eqAngle{i}(j,1) = min(thetaPlot{i}(:,j));
                    end
                elseif contains(robot.legObj{i}.joints{j+1},'TrF')
                    eqAngle{i}(j,1) = 0;
                elseif contains(robot.legObj{i}.joints{j+1},'CTr')
                    eqAngle{i}(j,1) = -35*pi/180;
                elseif contains(robot.legObj{i}.joints{j+1},'ThC2')
                    eqAngle{i}(j,1) = 0;
                elseif contains(robot.legObj{i}.joints{j+1},'ThC1')
                    eqAngle{i}(j,1) = 0;
                else %FTi
                    eqAngle{i}(j,1) = min(thetaPlot{i}(:,j))-(20*pi/180);
                end
            end
        end
        %%
        %Also figure out how many total equations we'll need to solve
        %everything
        matC = 0;
        matR = 0;
        for i=1:numLegs
            segM{i} = robot.legObj{i}.mass;

            matC = matC + 3*(robot.legObj{i}.numBodies-2+1 + robot.legObj{i}.numBodies-2);
            matR = matR + 6*(robot.legObj{i}.numBodies-2);
        end

        %Add in the equations for the thorax into the eq count
        matR = matR+6;

        for n=timeVector
            %First, set up the blank matrices we'll be populating for this
            %timestep
            A{n} = zeros(matR,matC);
            b{n} = zeros(matR,1);
            if n==1
                if orgLegMassTot < 1e-2
                    typX = 1e-6*ones(matC,1);

                else
                    typX = ones(matC,1);
                end
            end
            row = 1;
            legRowStart = 1;
            legColStart = 1;
            col = legColStart;
            jointLetters = {'F','A','B','C','D','E','H'};
            dirLetters = {'x', 'y', 'z'};
            swingCols{n} = [];
            for i=1:numLegs
                %Find the force and moment equations for each leg and assemble them into a
                %matrix
                segNum = robot.legObj{i}.numBodies-2; %Determine number of segments
                legEqNum = 6*segNum; %There will be 3 force equations and 3 moment equations for each segment
                legVarNum = 3*(segNum+1+segNum); %Each leg has segNum+1 forces and segNum moments (no moment at tarsus tip)
                torqueStartRow = legEqNum/2 + row; %The equations for moments will start halfway down the total equations for the leg
                torqueStartCol = legColStart + (segNum+1)*3; %Moment variables will also start after all of the force columns for the leg

                if n==1
                    numSegs(i) = segNum;
                    TvarStart{i} = torqueStartCol;
                    FvarStart{i} = legColStart;
                    FvarEnd{i} = torqueStartCol-1;
                    if i~=1
                        TvarEnd{i-1} = legColStart-1;
                    end
                    if i==6
                        TvarEnd{i} = matC;
                    end
                end

                for j=1:segNum %For each segment...
                    if j == 1 %If this is the tarsus...
                        if any(n == stanceIDs{i}) %Check if the leg is in stance
                            %If it is, populate the appropriate
                            %portions of the matrix w/ 1s
                            A{n}(row:row+2,col:col+2) = eye(3,3);
                            %Also find the r vector from the tarsus
                            %CoM to the tarsus tip
                            r_GLat = jointPoints{n}{i}(:,end-j+1)-comPoints{n}{i}(:,end-j);
                        else
                            swingCols{n} = [swingCols{n} col];
                            r_GLat = [0;0;0]; %If the leg is in swing, populate r with zeros
                        end
                        %Populate the r skew matrix for the foot
                        %force in the moment eqs
                        A{n}(torqueStartRow+3*(j-1):torqueStartRow+2+3*(j-1),(legColStart+3*(j-1)):(legColStart+2+3*(j-1))) = hat(r_GLat);
                        %Also populate the 1s for the TiTar moment
                        A{n}(torqueStartRow:torqueStartRow+2,torqueStartCol:torqueStartCol+2) = eye(3,3);
                    else %If this is not the tarsus...
                        A{n}(row+3*(j-1):row+3*(j-1)+2,col+3*(j-1):col+3*(j-1)+2) = -eye(3,3); %Populate the part of the equation for the lateral reaction forces
                        %Calculate the r vector from the body CoM
                        %to the lateral joint
                        r_GLat = jointPoints{n}{i}(:,end-j+1)-comPoints{n}{i}(:,end-j);
                        %Populate A with the r skew matrix for the
                        %lateral force
                        A{n}(torqueStartRow+3*(j-1):torqueStartRow+2+3*(j-1),(legColStart+3*(j-1)):(legColStart+2+3*(j-1))) = -hat(r_GLat);
                        %Populate for the reaction torque from the
                        %lateral joint
                        A{n}(torqueStartRow+3*(j-1):torqueStartRow+3*(j-1)+2,torqueStartCol+3*(j-2):torqueStartCol+3*(j-2)+2) = -eye(3,3);
                        %Populate the torques from the medial joint
                        A{n}(torqueStartRow+3*(j-1):torqueStartRow+3*(j-1)+2,torqueStartCol+3*(j-2)+3:torqueStartCol+3*(j-2)+3+2) = eye(3,3);
                    end
                    A{n}(row+3*(j-1):row+3*(j-1)+2,col+3+3*(j-1):col+3+3*(j-1)+2) = eye(3,3); %Populate the forces from the medial joint
                    %Calculate the r vector from the CoM to the
                    %medial joint
                    r_GMed = jointPoints{n}{i}(:,end-j) - comPoints{n}{i}(:,end-j);
                    %Populate the r skew matrix for the medial
                    %forces
                    A{n}(torqueStartRow+3*(j-1):torqueStartRow+2+3*(j-1),(legColStart+3*(j-1))+3:(legColStart+2+3*(j-1))+3) = hat(r_GMed);

                    if j == segNum %If we're on the last segment...
                        %Populate the appropriate cols at the very
                        %bottom of the equation for the thorax
                        A{n}(end-5:end-3,col+3+3*(j-1):col+3+3*(j-1)+2) = -eye(3,3); %ThC reaction forces
                        A{n}(end-2:end,(legColStart+3*(j-1))+3:(legColStart+2+3*(j-1))+3) = -hat(r_GMed); %Moments generated by reaction forces
                        A{n}(end-2:end,torqueStartCol+3*(j-2)+3:torqueStartCol+3*(j-2)+5) = -eye(3,3); %ThC reaction moments
                            
                    end

                    %------ Also populate the b matrix
                    %First, the ma portion of the force equations
                    b{n}(row+3*(j-1):row+3*(j-1)+2) = segM{i}(end-(j)) * segLinAcc{n}{i}(:,end-(j-1));
                    %For the y equation, add in the effect of gravity
                    b{n}(row+3*(j-1)+1) = b{n}(row+3*(j-1)+1) + segM{i}(end-(j)) * 9.81;
                    %Next, the moment equations
                    b{n}(torqueStartRow+3*(j-1):torqueStartRow+3*(j-1)+2) = cross(segAngVel{n}{i}(:,end-(j-1)),I_leg{i}{end-(j-1)}*segAngVel{n}{i}(:,end-(j-1))) ...
                        + I_leg{i}{end-(j-1)}*segAngAcc{n}{i}(:,end-(j-1)); %d/dt(H) term
                    if n == 1
                        if orgLegMassTot < 1e-2
                            typX(torqueStartRow+3*(j-1):torqueStartRow+3*(j-1)+2) = 10e-9*ones(3,1);
                        end
                    end

                    if robot.legObj{i}.motorEnabled(end-j) %If the joint has a servo actuating it along the omega axis...
                        b{n}(torqueStartRow+3*(j-1):torqueStartRow+3*(j-1)+2) = b{n}(torqueStartRow+3*(j-1):torqueStartRow+3*(j-1)+2) ...
                            + (springK * (theta{n}{i}(end-(j-1))-eqAngle{i}(end-(j-1))))*omegas{n}{i}(:,end-(j-1)); %Restoration torque from parallel springs
                    end
                    
                if i == numLegs && j == segNum %If we're on the last segment of the last leg
                    %Populate the last rows of the b matrix for the thorax
                    %First, the ma portion
                    b{n}(end-5:end-3) = segM{1}(1) * segLinAcc{n}{1}(:,1);
                    %For the y equation, add in gravity effects
                    b{n}(end-4) = b{n}(end-4) + segM{1}(1) * 9.81;
                    %Next, the moment equations
%                     b{n}(end-2:end) = cross([0;0;0],I_body*segAngVel{n}{1}(:,1)) ...
%                         + I_body*[0;0;0];
                    if n == 1
                        if orgLegMassTot < 1e-2
                            typX(end-5:end-3) = 10e-6*ones(3,1);
                        else
                            typX(end-5:end-3) = ones(3,1);
                        end
                    end
                end
                end
                %Once the leg is populated, define the next leg's
                %starting positions
                legRowStart = legRowStart + legEqNum;
                legColStart = legColStart+legVarNum;
                col = legColStart;
                row = legRowStart;
            end
            %Then add in some additional equations to explictly declare
            %that the F_F values in the swing legs are zero
            %First, check if it think all the legs are in swing
            if length(swingCols{n}) == 6
                for l = 1:length(swingCols{n-1})
                    B = swingCols{n} == swingCols{n-1}(l);
                    B = ~B;
                    swingCols{n} = B.*swingCols{n};
                end
                swingCols{n} = nonzeros(swingCols{n})';
            end
            for q = 1:length(swingCols{n})
                A{n} = [A{n}; zeros(3,matC)];
                A{n}(end-2:end,swingCols{n}(q):swingCols{n}(q)+2) = eye(3,3);
                b{n} = [b{n}; zeros(3,1)];
            end

        end

        %         %Create a symbolic x matrix to check if the equations look right
        %         row = 1;
        %         legRowStart = 1;
        %         legColStart = 1;
        %         col = legColStart;
        %         for i=1:numLegs
        %             segNum = robot.legObj{i}.numBodies-2;
        %             legVarNum = 3*(segNum+1+segNum); %Each leg has segNum+1 forces and segNum moments (no moment at tarsus tip)
        %             torqueStartRow = legEqNum/2 + row; %The equations for moments will start halfway down the total equations for the leg
        %             torqueStartCol = legColStart + (segNum+1)*3;
        %             for j=1:segNum+1
        %                 for n=1:3
        %                     x(col+3*(j-1)+(n-1),1) = sym(['F_' jointLetters{j} num2str(i) '_' dirLetters{n}]);
        %                 end
        %             end
        %             for j=1:segNum
        %                 for n=1:3
        %                     x(torqueStartCol+3*(j-1)+(n-1),1) = sym(['T_' jointLetters{j+1} num2str(i) '_' dirLetters{n}]);
        %                 end
        %             end
        %             legRowStart = legRowStart + legEqNum;
        %             legColStart = legColStart+legVarNum;
        %             col = legColStart;
        %             row = legRowStart;
        %         end

        %Then do the optimization to solve these massive matrices
        opts3 = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',20000,'UseParallel','always','ConstraintTolerance',1e-10,'StepTolerance',1e-10,'OptimalityTolerance',1e-10,'TypicalX',typX);
        for n=timeVector
            if n==1
                x0 = zeros(matC,1);
            else
                x0 = x{n-1};
            end
            if orgLegMassTot < 1e-2
                CEfun = @(x) 1e10*calculateCEeq(x,actuator,robot,omegas,n,TvarStart,TvarEnd,springK);
                x{n} = fmincon(CEfun,x0,[],[],1e10*A{n},1e10*b{n},[],[],[],opts3);
            else
                CEfun = @(x) calculateCEeq(x,actuator,robot,omegas,n,TvarStart,TvarEnd,springK);
                x{n} = fmincon(CEfun,x0,[],[],A{n},b{n},[],[],[],opts3);
            end
%             if n == 1
%             w = warning('query','last');
%             id = w.identifier;
%             warning('off',id)
%             end
            %Separate the forces and torques and figure out:
            % 1. How much torque is actually going to each actuator
            % 2. The components of each force in the leg's reference frame
            for i=1:robot.numLegs
                for j=1:numSegs(i)+1
                    FjointGlobal{i}{j}(:,n) = x{n}(FvarEnd{i}-(3*(j-1)+2):FvarEnd{i}-(3*(j-1)));
                    Fjoint{i}{j}(:,n) = g_b{n}(1:3,1:3,j+1)' * FjointGlobal{i}{j}(:,n);
                    if j ~= numSegs(i)+1
                        TjointGlobal{n}{i}(:,j) = x{n}(TvarEnd{i}-(3*(j-1)+2):TvarEnd{i}-(3*(j-1)));
                        Tjoint{i}{j}(:,n) = g_b{n}(1:3,1:3,j+1)' * TjointGlobal{n}{i}(:,j);
                        Tactuator{i}(n,j) = dot(omegas{n}{i}(:,j),TjointGlobal{n}{i}(:,j));
                    end
                end
            end
            fprintf([num2str(n), '.']);
        end
        fprintf('\n');

        timeVecPlot = linspace(0,stepPeriod,sampsPerStep);
        for i=1:numLegs
            %Plot the torques on each actuator
            tfig{i} = tiledlayout(2,1);
            title(tfig{i},{legNames{i}, springKPlotTitle, swingPeriodPlotTitle, configName, terrainShape},'Interpreter','tex')
            nexttile(1)
            for j=1:size(Tactuator{i},2)
                plot(timeVecPlot,Tactuator{i}(:,end-(j-1)))
                hold on
            end
%             if i > 4
%                 axis([0 timeVector(end) -2.7 1.3])
%             else
%                 axis([0 timeVector(end) -1.5 1.3])
%             end
%             yline(-1.2,'--')
%             yline(1.2,'--')
            ylabel('Torque (Nm)')
            xlim([0 stepPeriod])
            if contains(configName,'Drosophibot')
                yline(-1.5,'--')
                yline(1.5,'--')
            end
            nexttile(2)
            for j=1:size(Tactuator{i},2)
                plot(timeVecPlot,thetaPlot{i}(:,end-(j-1)))
                hold on
            end
            axis([0 stepPeriod -1.5 1.5])
            ylabel('Angle (rad)')
            if contains(configName, '2')
                if i < 5
                    legend('TiTar','FTi','TrF','CTr','ThC')
                else
                    legend('TiTar','FTi','TrF','CTr','ThC3','ThC1')
                end
            elseif contains(configName, 'Giga')
                legend('TiTar','FTi','TrF','CTr','ThC3', 'ThC1', 'ThC2')
            elseif contains(configName, 'Hex')
                legend('TiTar','FTi','CTr','ThC')
            else
                legend('TiTar','FTi','TrF','CTr','ThC3','ThC1','ThC2')
            end
            saveas(tfig{i},[savePathDir '\' legNames{i} '_Torques']);

            %Plot the GRF on each leg
            figure
            ffig{i} = tiledlayout(3,1);
            title(ffig{i},'Ground Reaction Forces, Global Frame','Interpreter','tex', 'FontWeight', 'bold')
            subtitle(ffig{i},{legNames{i}, [springKPlotTitle ', ' swingPeriodPlotTitle ', ' terrainShape], configName},'Interpreter','tex')
            nexttile(1)
            %Plot the x component of each of the forces
            plot(FjointGlobal{i}{end}(1,:))
            ylabel('X Component (N)')
            xlim([0 timeVector(end)])

            nexttile(2)
            plot(FjointGlobal{i}{end}(2,:))
            ylabel('Y Component (N)')
            xlim([0 timeVector(end)])
            
            nexttile(3)
            plot(FjointGlobal{i}{end}(3,:))
            ylabel('Z Component (N)')
            xlim([0 timeVector(end)])

            saveas(ffig{i},[savePathDir '\' legNames{i} '_GRFs']);

        end

        %Plotting for Brian
        if contains(configName, 'Male')
            if ~exist([savePathBase '\' saveName '\' springKPlot '\' swingPeriodPlot '\BrianData'])
                mkdir([savePathBase '\' saveName '\' springKPlot '\' swingPeriodPlot '\BrianData'])
            end
            savePathBrian = [savePathBase '\' saveName '\' springKPlot '\' swingPeriodPlot '\BrianData'];
            for i=1:numLegs
                for j=1:5
                    if j==1
                        FjointScale{i}{j} = Fjoint{i}{1} + Fjoint{i}{2} + Fjoint{i}{3};
                        TjointScale{i}{j} = Tjoint{i}{1} + Tjoint{i}{2} + Tjoint{i}{3};
                    else
                        FjointScale{i}{j} = Fjoint{i}{j+2};
                        TjointScale{i}{j} = Tjoint{i}{j+2};
                    end
                    for k=1:n
                        FjointScaleMags{i}{j}(k) = norm(FjointScale{i}{j}(:,k));
                        FjointScaleMags{i}{6}(k) = norm(Fjoint{i}{end}(:,k));
                        TjointScaleMags{i}{j}(k) = norm(TjointScale{i}{j}(:,k));
                    end
                    FjointScaleMax(i,j) = max(FjointScaleMags{i}{j});
                    TjointScaleMax(i,j) = max(TjointScaleMags{i}{j});
                end
                FjointScaleMax(i,6) = max(FjointScaleMags{i}{6});
            end
            jointNames = {'ThC','CTr','TrF','FTi','TiTar'};
            colors = {'#0072BD','#D95319','#EDB120'};
            timeVectorS = linspace(0,swingDuration/stanceDuty,n);
            for i=1:numLegs 
                for j=1:6
                    if j <= 5
                    figure
                    tl{i,j} = tiledlayout(3,2,'TileSpacing','compact','Padding','tight');
                    title(tl{i,j},[legNames{i} ' ' jointNames{j}],'Interpreter','tex','FontWeight', 'bold')
                    subtitle(tl{i,j},[springKPlotTitle ' ' '^{Nm}/_{rad}' ',   ' 'T_{step}= ' num2str(swingDuration/stanceDuty) ' s'],'Interpreter','tex')
                    for dir=1:3
                        nexttile
                        plot(timeVectorS,FjointScale{i}{j}(dir,:),'LineWidth',1.5,'Color',colors{dir})
                        grid on
                        xlim([0 timeVectorS(end)])
                        xlabel('Time (s)')
                        if dir==1
                            ylabel({'X Component (N)','Post <-- --> Ant'})
                            title('Forces')
                        elseif dir==2
                            ylabel({'Y (Lateral) Component (N)'})
                        else
                            ylabel({'Z Component (N)','Ventral <-- --> Dorsal'})
                        end

                        nexttile
                        plot(timeVectorS,TjointScale{i}{j}(dir,:),'LineWidth',1.5,'Color',colors{dir})
                        grid on
                        xlim([0 timeVectorS(end)])
                        xlabel('Time (s)')
                        if dir==1
                            ylabel({'X Component (Nm)','Post <-- --> Ant'})
                            title('Torques')
                        elseif dir==2
                            ylabel({'Y (Lateral) Component (Nm)'})
                        else
                            ylabel({'Z Component (Nm)','Ventral <-- --> Dorsal'})
                        end
                    end
                    saveas(tl{i,j},[savePathBrian '\' legNames{i} '_' jointNames{j}])
                    else
                        figure
                        tlGR{i} = tiledlayout(3,1);
                        title(tlGR{i},[legNames{i},' Tarsus Tip'],'Interpreter','tex', 'FontWeight', 'bold')
                        subtitle(tlGR{i},[springKPlotTitle ' ' '^{Nm}/_{rad}' ',   ' 'T_{step}= ' num2str(swingDuration/stanceDuty) ' s'],'Interpreter','tex')
                        for dir=1:3
                            nexttile
                            plot(timeVectorS,Fjoint{i}{end}(dir,:),'LineWidth',1.5,'Color',colors{dir})
                            xlim([0 timeVectorS(end)])
                            xlabel('Time (s)')
                            grid on
                            if dir==1
                                 ylabel({'X Component (N)','Post <-- --> Ant'})
                            elseif dir==2
                                 ylabel({'Y (Lateral) Component (N)'})
                            else
                                 ylabel({'Z Component (N)','Ventral <-- --> Dorsal'})
                            end
                        end
                        saveas(tlGR{i},[savePathBrian '\' legNames{i} '_Tarsus Tip'])
                    end

                end
            end
        end

        %Save the entire workspace for later
        save([savePath '\allVars'])

    end
end
end






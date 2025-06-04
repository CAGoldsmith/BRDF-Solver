classdef kinematicOrganism < matlab.mixin.SetGet
    %Kinetic_organism class to represent an organism with a number of legs
    %   Class has a cell array of legs, and a function for operating on all
    %   of the legs
    
    properties
        projFile = 'F:\My Documents\Animatlab\SingleJointImpedanceControl-PosStiffControl20140216-lengthmap\SingleJointImpedanceControl-PosStiffControl20140216-lengthmap.aproj';
        organismName = 'Organism_1';
        isRobot = true; %If isRobot, then muscles are ignored.
        updateFrequency;
        legObj = {}; % an array of the legs of the organism. Allows for multiple legs
        mirrored;
        numLegs;
        pJoints;
        toPlot;
        plastic;
        actuator;
        cpgProperties = [5,1,1,110,0,2,.05,20,1000,-0.05,0,0,0];
        gCpgToInter;
        swingDuration;
        maxStanceDuration;
        maxJointVelocity;
        maxDutyCycle;
        bodyTranslation;
        bodyRotation;
        constantAngleJoints;
        numTrans = 7;
        numRot = 15;
        ctr;
        delEdep = 220;
        delEhyp = -100;
        Rnet = 20;
    end
    
    methods
        function obj = kinematicOrganism(projFile, name, bodies, joints, ctr,...
                swingDuration, maxStanceDuration, bodyTranslation, bodyRotation,...
                constantAngleJoints, maxJointVelocity, mirrored, isRobot, toPlot, loopParameters, legPlastic, bodyPlastic, actuator)
            
            constructor = tic;
            
            obj.projFile = projFile;
            obj.organismName = name;
            
            obj.mirrored = mirrored;
            obj.isRobot = isRobot;
            obj.actuator = actuator;
            obj.plastic = bodyPlastic;
            
            obj.numLegs = size(bodies,2);
            obj.legObj = cell(obj.numLegs,1);
            obj.swingDuration = abs(swingDuration);
            obj.maxStanceDuration = abs(maxStanceDuration);
            obj.maxJointVelocity = abs(maxJointVelocity);
            if obj.maxStanceDuration < obj.swingDuration
                error('maxStanceDuration must be greater than swingDuration.')
            end
            obj.maxDutyCycle = obj.maxStanceDuration/(obj.maxStanceDuration+obj.swingDuration);
            obj.toPlot = toPlot;
            
            obj.bodyTranslation = abs(bodyTranslation);
            obj.bodyRotation = abs(bodyRotation);
            if size(constantAngleJoints,2) == 3
                obj.constantAngleJoints = constantAngleJoints;
            elseif isempty(constantAngleJoints)
                obj.constantAngleJoints = [-1,-1,-1];
            else
                error('constantAngleJoints must have exactly 3 columns.')
            end
            %First, make sure all "bodies" chains start with the same body,
            %which will be the main body segment.
            
            if obj.numLegs > 1
                for i=2:obj.numLegs
                    if ~strcmp(bodies{1,i-1},bodies{1,i})
                        error('All body chains must begin with the same root body.')
                    end
                end
            else
                
            end
            
            %ctr tells this function which joint is the primary "lifter"
            %joint in the leg. This joint's activity will be a function of
            %the remaining joint angles, to ensure a level stepping motion.
            if isscalar(ctr)
                obj.ctr = ctr + zeros(obj.numLegs,1);
            elseif length(ctr) == obj.numLegs
                obj.ctr = ctr;
            else
                error('ctr must be a scalar or a vector of length %i.',obj.numLegs)
            end
            
            for i=1:obj.numLegs
                %Each column represents another leg.
                try tempJoints = joints(:,i);
                catch err
                    errStr = strcat('The joint list has thrown the error ',err.identifier);
                    disp(errStr);
                end
                
                %Before this step, CHECK SIZES OF BODIES, JOINTS, MAX_SPEED
                if isempty(loopParameters)
                    obj.legObj{i} = leg(projFile,name,i,bodies(:,i),tempJoints,obj.ctr(i),obj.swingDuration,obj.maxStanceDuration,obj.maxDutyCycle,maxJointVelocity,isRobot,toPlot,legPlastic);
                else
                    obj.legObj{i} = leg(projFile,name,i,bodies(:,i),tempJoints,obj.ctr(i),obj.swingDuration,obj.maxStanceDuration,obj.maxDutyCycle,maxJointVelocity,isRobot,toPlot,legPlastic,loopParameters);
                end
                disp(['leg ',num2str(i),' constructed.']);
            end
            
            tConstruction = toc(constructor)/60;
            fprintf('It took %f minutes to construct the organism object.\n',tConstruction);
            
        end
        
        % add the leg object
        % returns the full list
        % NOTE use this, because the set_access is private for the property
        function [legList] = addLegObj( obj , legObj )
            obj.leg = [ obj.leg , legObj ];
            legList = obj.leg;
        end % addLeg w/ object
        
        function success = setLowLevelParams(obj)
            fprintf('Beginning to set organism parameters.\n')
            setParams = tic;
            
            fprintf('Beginning Jacobian computation.\n')
            obj.computeAllJacobians([],true);
%             fprintf('Beginning max torque computation.\n')
%             obj.computeAllMaxTorques();
            
%             if ~obj.isRobot
%                 fprintf('Beginning muscle force computation.\n')
%                 obj.computeAllMaxMuscForces();
%                 fprintf('Beginning joint control network training.\n')
%                 obj.trainAllJointControlNetworks();
%             else
% %                 
%             end
            
            paramTime = toc(setParams)/60;
            fprintf('It took %f minutes to set all the parameters for this organism.\n',paramTime);
            success = 1;
        end
        
        function [success,footPoints,omegas] = computeAllJacobians(obj,configs,combineFigures)
            
            if combineFigures && obj.toPlot
                if islogical(combineFigures)
                    kinematicFigure = figure;
                else
                    kinematicFigure = figure(combineFigures);
                end
            else
                kinematicFigure = [];
            end
            
            for i=1:obj.numLegs
                if iscell(configs)
                    [~,footPoints{i},obj.pJoints{i},omegas{i}] = obj.legObj{i}.computeJacobian(configs{i},kinematicFigure);
                elseif isempty(configs)
                    if isempty(obj.legObj{i}.restingConfig)
                        [~,footPoints{i},obj.pJoints{i},omegas{i}] = obj.legObj{i}.computeJacobian([],kinematicFigure);
                    else
                        [~,footPoints{i},obj.pJoints{i},omegas{i}] = obj.legObj{i}.computeJacobian(obj.legObj{i}.restingConfig,kinematicFigure);
                    end
                else
                    [~,footPoints{i},obj.pJoints{i},omegas{i}] = obj.legObj{i}.computeJacobian(configs(:,i),kinematicFigure);
                end
            end
            
            success = 1;
        end
        
        function success = computeAllMaxTorques(obj)
            allJoints = 1:obj.numLegs;
            if isempty(obj.mirrored)
                jointSet = allJoints;
            else
                jointSet = setdiff(allJoints,obj.mirrored(:,2));
            end
            
            for i=jointSet
                %Compute the parameters that can be mirrored for one
                %leg in each pair.
                obj.legObj{i}.computeMaxTorques(obj.bodyWeight);
                
                try indToMirror = obj.mirrored(obj.mirrored(:,1) == i,2);
                catch
                    indToMirror = 0;
                end
                
                %disp(['For leg ',num2str(i),' ind_to_mirror is ',num2str(ind_to_mirror)])
                if indToMirror
                    for j=2:obj.legObj{i}.numBodies
                        obj.legObj{indToMirror}.jointObj{j}.assignJointControllerProps(obj.legObj{i}.jointObj{j}.returnJointControllerProps());
                    end
                end
            end
            success = 1;
        end
        
        function success = computeAllMaxMuscForces(obj)
            allJoints = 1:obj.numLegs;
            if isempty(obj.mirrored)
                jointSet = allJoints;
            else
                jointSet = setdiff(allJoints,obj.mirrored(:,2));
            end
            
            for i=jointSet
                %Compute the parameters that can be mirrored for one
                %leg in each pair.
                obj.legObj{i}.computeMaxMuscleForce();
                
                try indToMirror = obj.mirrored(obj.mirrored(:,1) == i,2);
                catch
                    indToMirror = 0;
                end
                
                if indToMirror
                    for j=2:obj.legObj{i}.numBodies
                        obj.legObj{indToMirror}.jointObj{j}.assignJointControllerProps(obj.legObj{i}.jointObj{j}.returnJointControllerProps());
                    end
                end
            end
            success = 1;
        end
        
        function success = computeAllRobotAdapterGains(obj,jointHardwareLimits,savePath)
            saveLoc = [savePath '\Adapter Gains'];
            if ~exist(saveLoc)
                mkdir(saveLoc)
            end

            for i=1:obj.numLegs
                obj.legObj{i}.computeRobotAdapterGains(jointHardwareLimits,saveLoc);
            end
            success = 1;
        end
        
        function [configs,footPoints] = findRestingPosture(obj,jointVecs,standingHeight,boundingFactor,postureScaleFactor,terrainShape,ballParams)
            configs = cell(1,obj.numLegs);

            % FOR THIS CODE: Run the flat plane stuff first, then if necessary solve for
            % the ball position using flat plane as a starting position
            for i=1:obj.numLegs
                %define convenience handles for the legs and joints
                leg = obj.legObj{i};
                jointVec = jointVecs{i};

                %find the hardware limits of rotation, the range of motion,
                %and the number of joints.
                lb = leg.jointHardwareLimits(jointVec,1);
                ub = leg.jointHardwareLimits(jointVec,2);
                ran = ub - lb;
                n = length(ran);

                %find the mapping between x, the normalized variable, and
                %theta, the joint rotation. x = A*theta + b.
                %A = diag(1./ran);
                Ainv = diag(ran); %because the matrix is diagonal;
                b = -lb./ran;

                %the initial configuration is just 0.5, that is, halfway
                %through the range of motion.
                initConfig = 0.5+zeros(n,1);

                %find a scale factor that will ensure a parabola, g, that is in [-1,0] and is zero at x = -1 and 1.
                %this gives us a normalized objective function; inputs vary
                %Parabola in which g is 0 at 1 and -1
                scaleFac = 1/(n/4); %1/( (0.5+zeros(n,1))'*(0.5+zeros(n,1)) );
                gDistance = @(x) x'*(scaleFac*eye(n))*(x-ones(size(x))); %g = x(x-1)

                %We need an equality constraint that ensures the foot has
                %the proper height. We transform x into theta, and plug it
                %into the legHeightCon function.
                %legHeightCon is a separate file.
                fHeightCon = @(x) legHeightCon(Ainv*(x - b),leg,jointVec,standingHeight(i),postureScaleFactor,terrainShape);

                Aeqcon = [];
                beqcon = [];

                opts = optimoptions('fmincon', 'StepTolerance', 1e-10);
                %find the final normalized configuration.
                finalConfig = fmincon(gDistance,initConfig,[],[],Aeqcon,beqcon,boundingFactor+zeros(n,1),(1-boundingFactor)+zeros(n,1),fHeightCon,opts);
                %convert the normalized configuration to natural
                %coordinates.
                restingConfig = zeros(leg.numBodies-2,1);
                restingConfig(jointVec-1) = Ainv*(finalConfig - b);

                if contains(terrainShape, 'Ball')
                    [~,footPoint] = obj.legObj{i}.computeJacobian(restingConfig,[]);
                    newYCoord = (sqrt(ballParams.r^2 - (footPoint(1)-ballParams.center(1))^2 - (footPoint(3)-ballParams.center(3))^2) + ballParams.center(2));
                    footPointNew = [footPoint(1); newYCoord; footPoint(3)];
                    firstInd = jointVec(1)-1;
                    lastInd = jointVec(end)-1;
                    footCondition = @(vars) footCond(leg,vars,footPointNew,firstInd,lastInd);
                    angleEq = @(vars) norm(vars - restingConfig(firstInd:lastInd));
                    opts2 = optimoptions('fmincon', 'StepTolerance', 1e-10, 'ConstraintTolerance', 1e-10, 'Display','final');
                    ballRestConfig = fmincon(angleEq,restingConfig(firstInd:lastInd)+.1,[],[],Aeqcon,beqcon,[],[],footCondition,opts2);
                    restingConfig(jointVec-1) = ballRestConfig;
                else

                    %plot the normalized output relative to the bounds.
                    figure
                    plot(zeros(n,1),1:n,'r*')
                    hold on
                    plot(1+zeros(n,1),1:n,'g*')
                    plot(boundingFactor+zeros(n,1),1:n,'r--')
                    plot((1-boundingFactor)+zeros(n,1),1:n,'g--')
                    plot(finalConfig,1:n,'bo')
                    grid on
                    ylabel('joint number')
                    xlabel('rotation (normalized to range of motion)')
                    title(['Leg ',num2str(i)])
                end
                %set the leg's resting configuration.
                leg.setRestingConfig(restingConfig,terrainShape);
                configs{i} = restingConfig;
            end

            obj.setToPlot(true);
            [~,footPoints] = obj.computeAllJacobians(configs,true); %Automate here

        end
        
        function success = trainAllJointControlNetworks(obj)
            allJoints = 1:obj.numLegs;
            if isempty(obj.mirrored)
                jointSet = allJoints;
            else
                jointSet = setdiff(allJoints,obj.mirrored(:,2));
            end
            
            for i=jointSet
                %Compute the parameters that can be mirrored for one
                %leg in each pair.
                obj.legObj{i}.trainControlNetwork();
                
                try indToMirror = obj.mirrored(obj.mirrored(:,1) == i,2);
                catch
                    indToMirror = 0;
                end
                
                if indToMirror
                    for j=2:obj.legObj{i}.numBodies
                        obj.legObj{indToMirror}.jointObj{j}.assignJointControllerProps(obj.legObj{i}.jointObj{j}.returnJointControllerProps());
                    end
                end
            end
            success = 1;
        end
        
        function [rotationMat,translationMat,steppingKinematicsCell,stKin,swKin,pepJointAngles,aepJointAngles,pepConfigs,aepConfigs] = walkingInverseKinematics(obj,jointVecs,pBase,stanceDuty,stepHeight,stepProfile,terrainShape,ballParams)
            if contains(terrainShape, 'Ball')
                obj.numRot = 1;
                obj.numTrans = 2;
            end
            numTrials = obj.numRot*obj.numTrans;
            
            %Make sure we have a point on the robot's body to translate,
            %and pivot about.
            if isempty(pBase) && obj.numLegs >= 2
                %Define the position of a base point, the point between the
                %hind feet.
                pBase = (obj.pJoints{1}(:,end) + obj.pJoints{2}(:,end))/2;
            else
                if ~isequal(size(pBase),[3,1])
                    error('Base point for body transformation must be a 3 x 1 vector.')
                end
            end
            
            %Specify body translations and rotations with respect to the
            %base point.
            %Body translation deals with the spacial length of one step
            %Body rotation is how much the body can rotate with one step
            legTranslation = linspace(-obj.bodyTranslation,obj.bodyTranslation,obj.numTrans);
            if contains(terrainShape,'Ball')
                legRotation = 0;
            else
                legRotation = linspace(-obj.bodyRotation,obj.bodyRotation,obj.numRot);
            end
            [translationMat,~] = meshgrid(legTranslation,legRotation);
            for i = 1:obj.numTrans
                if legTranslation(i) < 0
                    thetaLimit = legTranslation(i) * obj.bodyRotation/obj.bodyTranslation + obj.bodyRotation;     
                else
                    thetaLimit = -legTranslation(i) * obj.bodyRotation/obj.bodyTranslation + obj.bodyRotation;
                end
                rotationMat(i,:) = linspace(-thetaLimit, thetaLimit,obj.numRot);
            end
            rotationMat = rotationMat';
            
            %Define numJoints to control looping for each leg's joints
            %synergies.
            if iscell(jointVecs) && (length(jointVecs) == obj.numLegs)
                numJoints = NaN(obj.numLegs,1);
                for i=1:obj.numLegs
                    numJoints(i) = length(jointVecs{i});
                end
            else
                error('jointVecs must be a cell array with one index for each leg.')
            end    
            %Using the ctr index, create a list of the "other" joints.
            numJoints = cellfun(@length,jointVecs)';
            othersVec = cell(obj.numLegs,1);
            for i=1:obj.numLegs
                othersVec{i} = jointVecs{i}; %1:numJoints(i);
                othersVec{i}(othersVec{i} == obj.ctr(i)) = [];
            end
            
            %Find the position of each foot wrt the base point. This will
            %be used to find the foot velocity for each trans and rot
            %velocity.
            rFootWrtBase = NaN(3,obj.numLegs,numTrials);
            
            %For each leg, for each translation, for each rotation, find
            %the velocity of the foot, and the joint synergy necessary to
            %move the foot in that direction.
            fff = figure;
            axis equal
            grid on
            xlabel('x')
            ylabel('z')
            drawnow
            
            %ffvecs = figure;
            
            %We will keep track of the joint angles for one complete
            %steppig cycle for each leg, for each body translation, for
            %each body rotation.
            steppingKinematicsCell = cell(obj.numLegs,obj.numTrans,obj.numRot);
            
%             %Joint angles as reconstructed by the controller. These will be
%             %used for comparison to the exact mathematical scenario.
%             reconstructedKinematicsCell = cell(obj.numLegs,obj.numTrans,obj.numRot);
            
%             %Reconstructed angles in which the CTr angle is a function of
%             %the other leg angle joints, based on the interjoint reflexes
%             %observed by Hess et al. 1997 and 1999.
%             reconstructedKinematicsHessCell = cell(obj.numLegs,obj.numTrans,obj.numRot);
            
            %Joint angles at the PEP and AEP of each step for each leg, for
            %each body translation, for each body rotation.
            pepJointAngles = cell(obj.numLegs,obj.numTrans,obj.numRot);
            aepJointAngles = cell(obj.numLegs,obj.numTrans,obj.numRot);

            pepConfigs = cell(obj.numLegs,1);
            aepConfigs = cell(obj.numLegs,1);

            kinDS = cell(obj.numLegs,obj.numTrans,obj.numRot);
            
            planarFitCoeffs = cell(obj.numLegs,1);
            
            numSampsPerStep = 41; %Must be an odd number.
            for i=1:obj.numLegs
                %define stanceDirections. This is one of our outputs. It
                %isn't used for computation, so its clunky form (cell of
                %matrices) is probably ok.
                
                %convenience handle
                jointVec = jointVecs{i};
                
                pepConfigs{i} = NaN(obj.numTrans,obj.numRot,obj.legObj{i}.numBodies-2);
                aepConfigs{i} = NaN(obj.numTrans,obj.numRot,obj.legObj{i}.numBodies-2);
                                
                %if this leg has a body joint, then flex it proportionally
                %to the heading angle
                for j=1:obj.numRot
                    %                         rFootWrtBase(:,i,j) = obj.pJoints{i}(:,end) - pBase;
                    [~,pFoot] = obj.legObj{i}.computeJacobian(obj.legObj{i}.getRestingConfig,[]);
                    rFootWrtBase(:,i,j) = pFoot - pBase;
                end
                
                for j=1:obj.numTrans
                    for k=1:obj.numRot
                        j_p = squeeze(rFootWrtBase(:,i,k)); %in spatial frame
                        r = [translationMat(k,j);0;0]; %in body frame; how far body moves with each step %CHECK HERE FIRST 
                        psi = rotationMat(k,j);
                        
                        R_oneFrame = [  cos(psi/(numSampsPerStep-1))     0   sin(psi/(numSampsPerStep-1));...
                                        0                   1   0;...
                                        -sin(psi/(numSampsPerStep-1))    0   cos(psi/(numSampsPerStep-1))];
                        R_oneFrameNeg = [ cos(-psi/(numSampsPerStep-1))    0   sin(-psi/(numSampsPerStep-1));...
                                    0                   1   0;...
                                    -sin(-psi/(numSampsPerStep-1))   0   cos(-psi/(numSampsPerStep-1))];
                        r_oneFrame = r/(numSampsPerStep-1);
                                
                        mid = (numSampsPerStep+1)/2;
                        qVec = NaN(3,mid); %qVec is position of the foot
                        qVec(:,mid) = j_p;
                        for m=mid+1:numSampsPerStep
                            qVec(:,m) = R_oneFrameNeg*qVec(:,m-1)-r_oneFrame;
                        end
                        for m=mid-1:-1:1
                            qVec(:,m) = R_oneFrame*(qVec(:,m+1)+r_oneFrame);
                        end
                        
                        qVec = qVec - j_p;
                        
                        if legTranslation(j) < 0 && legRotation(k) < 0
                            color = 'r';
                        elseif legTranslation(j) > 0 && legRotation(k) < 0
                            color = 'k';
                        elseif legTranslation(j) < 0 && legRotation(k) > 0
                            color = 'g';
                        elseif legTranslation(j) > 0 && legRotation(k) > 0
                            color = 'b';
                        else
                            color = 'm';
                        end
                        
                        figure(fff)
                        hold on
                        plot(j_p(1)+qVec(1,:),j_p(3)+qVec(3,:),color)
                        xlabel('x -> forward')
                        ylabel('z -> right')
%                         xlim([-.3,.3])
%                         ylim([-.3,.3])
                        axis equal
                        hold off
                                                
                        startConfig = obj.legObj{i}.getRestingConfig;
                        
                        stepFrequency = 1/(2*obj.swingDuration); %This could be changed to facilitate tetrapod, etc.
                        ThCPos = obj.pJoints{i}(:,1);

                        if contains(terrainShape, 'Ball')
                            %Modify qVec's y coordinate to be the right
                            %height for the ball
                            qVecGlobal = qVec + pFoot;
                            qVecGlobalBall = qVecGlobal;
                            for n=1:length(qVec)
                                qVecGlobalBall(2,n) = (sqrt(ballParams.r^2 - (qVecGlobal(1,n)-ballParams.center(1))^2 - (qVecGlobal(3,n)-ballParams.center(3))^2) + ballParams.center(2));
                            end
                            qVec = qVecGlobalBall - pFoot;
                        end
                        [steppingKinematicsCell{i,j,k},stKin{i,j,k},swKin{i,j,k},AEPID,PEPID,footPath] = obj.legObj{i}.steppingKinematics(startConfig,jointVec,qVec,numSampsPerStep,stepHeight(i),stanceDuty,stepFrequency,stepProfile,terrainShape);
                        pepJointAngles{i,j,k} = steppingKinematicsCell{i,j,k}(:,PEPID);
                        aepJointAngles{i,j,k} = steppingKinematicsCell{i,j,k}(:,AEPID);
                        pepConfigs{i}(j,k,:) = pepJointAngles{i,j,k} - obj.legObj{i}.restingConfig;
                        aepConfigs{i}(j,k,:) = aepJointAngles{i,j,k} - obj.legObj{i}.restingConfig;
                    end
                end
            end

            
            %max forward, max forward and turn, no forward and max turn
            fwd = [obj.numTrans,obj.numTrans,ceil(obj.numTrans/2)];
            trn = [ceil(obj.numRot/2),obj.numRot,obj.numRot];
            
            toAnimate = false;
            
%             for i=[1,4,6]
%                 if ~toAnimate
%                     footpathDiagram = false;
%                 else
%                     %What value do you want it to be?
%                     footpathDiagram = true;
%                 end
%                 
%                 obj.animateFullKinematics(steppingKinematicsCell(:,fwd(i),trn(i)),pBase,stepHeight,toAnimate,footpathDiagram,['Just IK ',num2str(i)]);
% 
%                 footpathDiagram = false;
% %                 groundPenRecon = obj.animateFullKinematics(reconstructedKinematicsCell(:,fwd(i),trn(i)),pBase,stepHeight,toAnimate,footpathDiagram,['Linear Interp of AEP to PEP ',num2str(i)]);
% %                 
% %                 groundPenReconHess = obj.animateFullKinematics(reconstructedKinematicsHessCell(:,fwd(i),trn(i)),pBase,stepHeight,toAnimate,footpathDiagram,['Linear Interp plus reflexes to CTr',num2str(i)]);
%                 
%             end
        end
        
        function groundPenetration = animateFullKinematics(obj,steppingKinematics,pBase,stepHeight,toAnimate,toPlotHistory,animationTitle)
            m = size(steppingKinematics{1},2);
            offset = zeros(obj.numLegs,1);
            offset([2,3,6]) = ceil(m/2);
            
            dsFactor = 5;
            numCycles = 3;
            
            frameVec = 1:dsFactor:numCycles*m;
            oneCycleFrameVec = 1:dsFactor:m;
            
            dt = 2*obj.swingDuration/length(oneCycleFrameVec);
            
            %should be n_legs, n_samples
            groundPenetration = NaN(length(steppingKinematics),length(oneCycleFrameVec));
                 
            xll = -.45;
            xul = .25;
%             yll = 0;
%             yul = .1;
%             zll = -.3;
%             zul = .3;
            yll = -.35;
            yul = .35;
            zll = -.1;
            zul = .3;
            
            if toAnimate
                hAnimate = figure;
                hTitle = title(animationTitle);
                view(140,40)
                axis equal
                hold on
                grid on
                xlims = [xll,xul];
                xlim(xlims)
                ylims = [yll,yul];
                ylim(ylims)
                zlims = [zll,zul];
                zlim(zlims)
                
                myVideo = VideoWriter(animationTitle);
                myVideo.FrameRate = 20;
                myVideo.Quality = 100;
                open(myVideo);
            end
            
            if toPlotHistory
                hStepHistory = figure;
                title(animationTitle);
                hold on
                view(-316,36)
                axis equal
                grid on
                colors = get(gca,'colororder');
                footPaths = cell(6,1);
                
                xlims = [xll,xul];
                xlim(xlims)
                ylims = [yll,yul];
                ylim(ylims)
                zlims = [zll,zul];
                zlim(zlims)
            end
            
            hJoints = cell(obj.numLegs,1);
            hFeet = cell(obj.numLegs,1);
            
            initialStep = true;
            toContinue = true;
            while toContinue
                k = 0;
                for i=frameVec
                    k = k + 1;
                    
                    hTitle.String = sprintf('%s t = %0.3f s',animationTitle,k*dt);
                    for j=1:obj.numLegs
                        [~,f_p,j_p] = obj.legObj{j}.computeJacobian(steppingKinematics{j}(:,mod(i + offset(j),m)+1),[]); %FIX ME
                            
                        if toAnimate
                            figure(hAnimate)
                            if initialStep
                                hJoints{j} = plot3(j_p(1,:),-j_p(3,:),j_p(2,:),'linewidth',1);
                                hFeet{j} = plot3(f_p(1),-f_p(3),f_p(2),'ko');
                                
                                [X,Y] = meshgrid(2*xlims,2*ylims);
                                Z = -stepHeight+zeros(size(X));
                                floorSurf = surf(X,Y,Z,'edgealpha',0,'facealpha',.2);
                                drawnow
                                
                                xlim(xlims);
                                ylim(ylims);
                            else
                                hJoints{j}.XData = j_p(1,:);
                                hJoints{j}.YData = -j_p(3,:);
                                hJoints{j}.ZData = j_p(2,:);

                                hFeet{j}.XData = f_p(1);
                                hFeet{j}.YData = -f_p(3);
                                hFeet{j}.ZData = f_p(2);
                            end
                        end
                        
                        if i < m
                            groundPenetration(j,k) = max(0,pBase(2) - f_p(2));
                        end
                         
                        if toAnimate
                            figure(hAnimate)
                            if f_p(2) > pBase(2) + 0.01*stepHeight
                                %in swing phase
                                hFeet{j}.Visible = 'off';
                            else
                                hFeet{j}.Visible = 'on';
                            end
                            
%                             limsChanged = false;
%                             if any(j_p(1,:) > xul)
%                                 xul = max(j_p(1,:));
%                                 limsChanged = true;
%                             end
%                             if any(j_p(1,:) < xll)
%                                 xll = min(j_p(1,:));
%                                 limsChanged = true;
%                             end
%                             if any(-j_p(3,:) > yul)
%                                 yul = max(j_p(2,:));
%                                 limsChanged = true;
%                             end
%                             if any(-j_p(3,:) < yll)
%                                 yll = min(j_p(2,:));
%                                 limsChanged = true;
%                             end
%                             if any(j_p(2,:) > zul)
%                                 zul = max(j_p(3,:));
%                                 limsChanged = true;
%                             end
%                             if any(j_p(2,:) < zll)
%                                 zll = min(j_p(3,:));
%                                 limsChanged = true;
%                             end
                            limsChanged = true;
                            if limsChanged
                                xlim([xll,xul])
                                ylim([yll,yul])
%                                 zlim([zll,zul])
                            end
                        end
                        
                        if toPlotHistory
    %                         if (i <= m/2+1) && mod(i-1,dsFactor*8) == 0
                            if (i <= m+1) && mod(i-1,dsFactor*8) == 0
                                if strcmp(hFeet{j}.Visible,'on')
                                    figure(hStepHistory)
                                    hold on
    %                                 pTemp = plot3(j_p(1,:),j_p(2,:),j_p(3,:),'linewidth',2,'color',colors(j,:));
                                    pTemp = plot3(j_p(1,:),-j_p(3,:),j_p(2,:),'linewidth',1,'color',colors(j,:));
                                    pTemp.Color(4) = .25 + .75*mod(i-1,m/2)/(m/2);
                                    footPaths{j} = [footPaths{j},[j_p(1,end);-j_p(3,end);j_p(2,end)]];
                                    drawnow
                                end
                            end 
                        end
                    end
                    
                    if initialStep
                        initialStep = false;
                    end
                    
                    if toAnimate
                        drawnow
                        figure(hAnimate);
                        frame = getframe(gcf); %get frame
                        writeVideo(myVideo, frame);
                    end
                end
                
                if toPlotHistory
                    %add ground plane
                    figure(hStepHistory)
                    ax = gca;
                    xl = ax.XLim;
                    yl = ax.YLim;
                    [X,Y] = meshgrid(xl,yl);
                    Z = -stepHeight+zeros(size(X));
                    surf(X,Y,Z,'edgealpha',0,'facecolor',.9+zeros(1,3));

                    for j=1:obj.numLegs
                        if isequal(footPaths{j}(:,1),footPaths{j}(:,end))
                            footPaths{j}(:,end) = [];
                        end
                        plot3(footPaths{j}(1,:),footPaths{j}(2,:),footPaths{j}(3,:),'k','linewidth',1)
                    end
                    
                    hStepHistory.Position = [865.8000  253.0000  286.7200  215.0400];
                    xticks([])
                    yticks([])
                    zticks([])                   
                    
                end
                    
                if toAnimate
                    close(myVideo);
                    
%                     response = input('Press [Enter] to repeat animation, enter any other character to continue. ','s');
                    response = 'a';
                    if strcmp(response,'')
                        %do nothing
                    else
                        toContinue = false;
                    end
                else
                    toContinue = false;
                end
            end
        end
        
        function out = sampleConfigs(obj,allConfigs,numTrials)
            wholeBodyConfig = cell(obj.numLegs,1);
            for i=1:10
                trial = randi(numTrials);
                for j=1:obj.numLegs
                    wholeBodyConfig{j} = allConfigs{j}(trial,:)';
                end
                figure
                obj.computeAllJacobians(wholeBodyConfig,true);
            end
            out = 123;
        end
        
        function [G, Vss] = computeCpgToMnSynapse(obj)
            addpath('CPG\CPGDESIGN')
            Erest = obj.cpgProperties(5);
            R = 20; %operating range for nonspiking signals in the network
            
            %First, obtain a CPG activation of -40 mV by choosing Sm and Sh
            fNprops = @(S) [obj.cpgProperties(1:6),S,obj.cpgProperties(8:9),-S,obj.cpgProperties(11:13)];
            vTonicDesired = Erest+R; %mV
            f = @(S)dV_dt( vTonicDesired,hinf_of_v(vTonicDesired,fNprops(S)),0,0,0,fNprops(S));
            Sm = fzero(f,[0,1]);
            Sh = -Sm;
            obj.cpgProperties(7) = Sm;
            obj.cpgProperties(10) = Sh;
            
            %Next, verify this result.
            f = @(x) dV_dt(x,hinf_of_v(x,obj.cpgProperties),0,0,0,obj.cpgProperties);
            Vss = fzero(f,[-20,60]);
            if abs( (Vss - vTonicDesired)/(Vss) ) > 0.10
                warning('CPG steady state voltage is more than 10\% wrong. This may affect accuracy.')
            end
            
            %Finally, find the synapse conductivity that will make the
            %postsynaptic output R when the presynaptic input is R.
            synFac = (Vss - Erest)/R;
            f = @(x) dV_dt(Erest+R,0,synFac*x,300,0,obj.cpgProperties);
            G = fzero(f,[0,1]);
            
            %Write these
            obj.gCpgToInter = G;
            for i=1:obj.numLegs
                for j=2:obj.legObj{i}.numBodies
                    obj.legObj{i}.jointObj{j}.setCpgParams(G,Vss);
                end
            end
            
        end
        
        function newPlotVal = setToPlot(obj, toPlot)
            obj.toPlot = toPlot;
            for i=1:obj.numLegs
                obj.legObj{i}.setToPlot(toPlot);
            end
            newPlotVal = obj.toPlot;
        end
        
        function success = forAll( obj , functionNameString , varargin )
            % give it a function to run, and it will run for all of the
            % legs in the object
            
            % passes the rest of the arguments as the cell array varargin
            % NOTE this might obscure some of the argument errors, I'm not
            % sure how carefully Matlab checks.
            success = 1;
            for i = 1:len(obj.leg)
                try
                    obj.leg(i).feval(functionNameString, varargin);
                catch exception
                    disp(['exception in function evaluation: ',exception]);
                    success = 0;
                end % end function callign try/catch sequence
            end % loop through all of the legs
        end % forAll function
        
        function newAprojFile = applyAllParamValues(obj)
            
            %If the user is modifying a "..._mod.aproj" type file, the
            %standalone file will not include the _mod (or any multiplicity
            %of them). Therefore we need to modify the project file name
            %(or a copy of it) and remove the "_mod"s to compute the actual
            %standalone file name.
            %             pseudo_proj_file = obj.projFile;
            %             while strcmp(pseudo_proj_file(end-9:end-6),'_mod')
            %                 pseudo_proj_file = [pseudo_proj_file(1:end-10),pseudo_proj_file(end-5:end)];
            %             end
            %
            %             %Here is the .asim file to find.
            %             asim_file = [pseudo_proj_file(1:end-6),'_Standalone.asim'];
            
            allJoints = 1:obj.numLegs;
            if isempty(obj.mirrored)
                jointSet = allJoints;
            else
                jointSet = setdiff(allJoints,obj.mirrored(:,2));
            end
            
            props = {};
            for i=1:obj.numLegs
                props = [props;obj.legObj{i}.listAllLegParams(any(jointSet == i))];
            end
            newAprojFile = obj.updateFileInfo(props,'AnimatlabProperties.mat');
            %             new_aproj_file = updateFileInfo(props,'AnimatlabProperties.mat',obj.legObj{1}.original_text);
            
        end
        
        function writtenFile = updateFileInfo( obj, infoCellArray, propertyFile )
            %UPDATE_FILE_INFO Last updated 17 Aug 16
            %   input:
            %   cell array
            %   object names | property name | numerical
            %
            %   output:
            %   take existing file and write new properties, return the file
            
            cellAprojText = obj.legObj{1}.originalText;
            numRows = size(infoCellArray,1);
            i = 1;
            
            while i <= numRows
                %                 if isequal(infoCellArray{i,1},'LM_ThC1 Trans') || isequal(infoCellArray{i,1},'LM_ThC1 Rotate')
                %                     obj.findID(infoCellArray{i,1},cellAprojText)
                %                     keyboard
                %                 end
                %         disp(['processing row ' num2str(row)]);
                % get the right id for the given object name
                [propToFind, specInfo, pageNumber] = obj.identifyProperty(propertyFile, infoCellArray{i,2});
                if ~strcmp(propToFind , 'Error')
                    specID = obj.findID(infoCellArray{i,1},pageNumber,cellAprojText);
                    if ~strcmp(specID, 'Error')
                        %with valid id + property name, process the value (check range
                        % + scale)
                        [val, sca, act, error] = obj.processNumerical(infoCellArray{i,3}, propToFind);
                        if ~error
                            % if the numerical processing doesn't throw an error:
                            % make the changes to the file for item id, property,
                            % value, scale, actual
                            if specInfo == 1
                                disp(['[object] ' infoCellArray{i,1} ...
                                    ' [id] ' specID ...
                                    ' [property] m' propToFind ...
                                    ' [val/sca/act] ' num2str(val) '/' sca '/' num2str(act)]);
                            elseif specInfo == 2
                                disp(['[object] ' infoCellArray{i,1} ...
                                    ' [id] ' specID ...
                                    ' [property] h' propToFind ...
                                    ' [val/sca/act] ' num2str(val) '/' sca '/' num2str(act)]);
                            elseif specInfo == 3
                                disp('damping')
                            else
                                disp(['[object] ' infoCellArray{i,1} ...
                                    ' [id] ' specID ...
                                    ' [property] ' propToFind ...
                                    ' [val/sca/act] ' num2str(val) '/' sca '/' num2str(act)]);
                            end
                            
                            if isempty(cellAprojText)
                                wasEmpty = true;
                            else
                                wasEmpty = false;
                            end
                            cellAprojTextNew = obj.setInfo(specID, propToFind, specInfo, val, sca, act, cellAprojText);
                            if (isempty(cellAprojTextNew) && ~wasEmpty)
                                %do nothing
                                warning(['The item "',infoCellArray{i,1},'" may not exist within the simulation. Please add it to the simulation, save it, and export another standalone simulation.']);
                            else
                                cellAprojText = cellAprojTextNew;
                            end
                            
                        end
                    else
                        disp(['[ERROR] No ID found for [object] ' infoCellArray{i,1}]);
                    end
                else
                    disp(['[ERROR] Invalid property for [object] ' infoCellArray{i,1} ...
                        ' [property] ' infoCellArray{i,2}]);
                end
                i = i + 1;
                
            end
            % write the new file to disk with the old data + changes
            disp('ran to file write');
            fid = fopen([obj.projFile(1:length(obj.projFile)-6) '_mod.aproj'], 'wt');
            fprintf(fid, '%s\n' , cellAprojText{:});
            fclose(fid);
            
            writtenFile = [obj.projFile(1:length(obj.projFile)-6) '_mod.aproj'];
        end

        function success = setMotionAngles(obj,jointVecs,pepConfigs,aepConfigs,rotationMat,translationMat,savePath)
            saveLoc1 = [savePath '\Body Motion Maps'];
            saveLoc2 = [savePath '\Synapse Strengths'];
            if ~exist(saveLoc1)
                mkdir(saveLoc1)
            end
            if ~exist(saveLoc2)
                mkdir(saveLoc2)
            end
            numJoints = NaN(obj.numLegs,1);
            synTableColNames = {'Rotation (nS)', 'Translation (nS)','PEP+ Rest Potential (mV)','PEP- Rest Potential (mV)'};
            for i=1:obj.numLegs
                numJoints(i) = length(jointVecs{i});
                synTableRowNames{i} = [];
                gScaleRot = zeros(numJoints(i),1);
                gScaleTrans = zeros(numJoints(i),1);
                aepRestingPot = zeros(numJoints(i),1);
                pepRestingPot = zeros(numJoints(i),1);
                for m=1:numJoints(i)
                    obj.legObj{i}.jointObj{jointVecs{i}(m)}.bodyTranslation = translationMat;
                    obj.legObj{i}.jointObj{jointVecs{i}(m)}.bodyRotation = rotationMat;
                    obj.legObj{i}.jointObj{jointVecs{i}(m)}.setLocomotionAngles(pepConfigs{i}(:,:,jointVecs{i}(m)-1),aepConfigs{i}(:,:,jointVecs{i}(m)-1));
                    obj.legObj{i}.jointObj{jointVecs{i}(m)}.designBodyMotionMap(saveLoc1);
                    gScaleRot(m) = obj.legObj{i}.jointObj{jointVecs{i}(m)}.gScaleRot;
                    gScaleTrans(m) = obj.legObj{i}.jointObj{jointVecs{i}(m)}.gScaleTrans;
                    aepRestingPot(m) = obj.legObj{i}.jointObj{jointVecs{i}(m)}.aepRestingPot;
                    pepRestingPot(m) = obj.legObj{i}.jointObj{jointVecs{i}(m)}.pepRestingPot;
                    synTableRowNames{i} = [synTableRowNames{i}, convertCharsToStrings(obj.legObj{i}.jointObj{jointVecs{i}(m)}.jointName)];
                end
                synLegTable{i} = table(gScaleRot*1e3,gScaleTrans*1e3,pepRestingPot,aepRestingPot,'VariableNames',synTableColNames,'RowNames',synTableRowNames{i});
                writetable(synLegTable{i},[saveLoc2,'\Leg ' num2str(i) '.txt'],'Delimiter',' ','WriteRowNames',1)
            end
            success = 1;
        end
    end
    
    methods(Static)
        
        function [fullID] = findID(objectName, pageNumber, cellAprojFile)
            fullID = 'Error';
            % while there are lines to lookup...
            i = 1;
            numRows = size(cellAprojFile,1);
            foundFlag = 0;
            
            while ~foundFlag && i <= numRows
                % if the current line is long enough to have a name...
                currentLine = cellAprojFile{i};
                if (length(currentLine)> length('<Name></Name>') && strncmpi(currentLine,'<Name>',6)) ||...
                        (length(currentLine)> length('<Text></Text>') && strncmpi(currentLine,'<Text>',6))
                    shortName = currentLine(7:length(currentLine)-7);
                    shortName2 = 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx';
                    if (length(currentLine)> length('<ModuleName></ModuleName>'))
                        shortName2 = currentLine(13:length(currentLine)-13);
                    end
                    % if the name is correct...
                    if strcmp(objectName,shortName) || strcmp(objectName,shortName2)
                        % now to find the id
                        %process backwards to the opening of the object tag
                        % (open tag w/o data)
                        objectStartFlag = 0;
                        while(~objectStartFlag && i>0)
                            currentLine = cellAprojFile{i};
                            numOpenTags = length(strfind(currentLine,'<'));
                            numCloseTags = length(strfind(currentLine,'>'));
                            numSlash = length(strfind(currentLine,'/'));
                            if (length(currentLine)> length('<>') && numOpenTags == 1 ...
                                    && numCloseTags == 1  && numSlash == 0)
                                % if the current line has just a tag
                                objectStartFlag = 1;
                            else
                                i = i -1;
                            end
                        end
                        % then move down or up until an <ID> is found.
                        %For synapses, the ID is at the bottom, so we need
                        %to move up.
                        if pageNumber == 2
                            iChange = -1;
                        else
                            iChange = 1;
                        end
                        
                        linksFlag = 0;
                        while ~foundFlag && i<length(cellAprojFile) && length(currentLine)>4 && (linksFlag || ~strncmpi(currentLine,'<ID>',4))
                            currentLine = cellAprojFile{i};
                            if strncmpi(currentLine,'<InLinks>',length('<InLinks>')) || strncmpi(currentLine,'<OutLinks>',length('<OutLinks>'))
                                linksFlag = 1;
                            elseif strncmpi(currentLine,'</InLinks>',length('</InLinks>')) || strncmpi(currentLine,'</OutLinks>',length('</OutLinks>'))
                                linksFlag = 0;
                            else
                                % otherwise, we have a valid id
                                if length(currentLine) > length('<ID></ID>') && strncmpi(currentLine,'<ID>',4)
                                    idLine = currentLine;
                                    fullID = idLine(5:length(idLine)-5);
                                    foundFlag = 1;
                                end
                            end
                            
                            i = i + iChange;
                        end
                    end
                end
                i = i+1;
            end
        end
        
        function [value, scale, actual, err] = processNumerical(actual, property)
            if isempty(actual)
                error([property,' may not be left blank.']);
            end
            
            actualCheck = str2double(actual);
            if isnan(actualCheck)
                value = actual;
                scale = [];
                err = 0;
            else
                actual = actualCheck;
                err = 0;
                
                if actual == 0
                    value = num2str(0);
                    scale = 'None';
                else
                    minorCount = 0;
                    tempAct = actual;
                    if abs(tempAct) > 1000
                        while abs(tempAct)>1000
                            tempAct = tempAct/1000;
                            minorCount = minorCount + 1;
                        end
                        value = num2str(tempAct);
                        if minorCount == 0
                            scale = 'None';
                        elseif minorCount == 1
                            scale = 'kilo';
                        elseif minorCount == 2
                            scale = 'mega';
                        elseif minorCount == 3
                            scale = 'giga';
                        elseif minorCount == 4
                            scale = 'tera';
                        else
                            scale = 'Error';
                        end
                    elseif abs(tempAct) <.1
                        while abs(tempAct)<.1
                            tempAct = tempAct*1000;
                            minorCount = minorCount + 1;
                        end
                        value = num2str(tempAct);
                        if minorCount == 0
                            scale = 'None';
                        elseif minorCount == 1
                            scale = 'milli';
                        elseif minorCount == 2
                            scale = 'micro';
                        elseif minorCount == 3
                            scale = 'nano';
                        elseif minorCount == 4
                            scale = 'pico';
                        else
                            scale = 'Error';
                        end
                    else
                        scale = 'None';
                        value = actual;
                    end
                end
            end
        end
        
        function [propToFind,specialInfo,pageNumber] = identifyProperty(propFile, lookupString)
            if strncmpi(lookupString,'m',1)
                specialInfo = 1;
            elseif strncmpi(lookupString,'h',1)
                specialInfo = 2;
            elseif strncmpi(lookupString,'damping',7)
                specialInfo = 3;
            else
                specialInfo = 0;
            end
            availableProperties = load(propFile);
            propSet = availableProperties.Properties;
            propToFind = 'Error';
            % lookup the property in the animatlab file
            numPages = length(availableProperties.Properties);
            for pageNumber = 1:numPages
                currentPage = propSet{pageNumber,1};
                sizeInfo = size(currentPage);
                sizeInfo = sizeInfo(1);
                for j = 1:sizeInfo % the length of the cell
                    if strcmp(currentPage{j,1},lookupString)
                        propToFind = currentPage{j,4};
                        return;
                    end
                end
            end
        end
        
        function fileWithEdit = setInfo(fullID, propertyName, modifier, value, scale, actual, celledProjFile)
            % Test
            % find the correct line by linear scan
            i = 1;
            fileWithEdit = [];
            replaced = 0;
            numRows = size(celledProjFile,1);
            
            mModFlag = 0;
            hModFlag = 0;
            skipEntry = 0;
            if modifier == 1
                mModFlag = 1;
            elseif modifier == 2
                hModFlag = 1;
            elseif modifier == 3
                skipEntry = 1;
            end
            while(i < numRows && ~replaced)
                if ischar(celledProjFile{i}) && length(celledProjFile{i}) > length('<ID></ID>')...
                        && strncmpi(celledProjFile{i},'<ID>',4)
                    
                    %identify the ID, if this line has one.
                    toCompare = celledProjFile{i};
                    toCompare = toCompare(5:length(toCompare)-5);
                    % if the correct ID is found
                    if strcmp(toCompare,fullID)
                        linksFlag = 0;
                        caA_m_flag = 0;
                        caD_h_flag = 0;
                        i = i+1;
                        % cycle through the properties until it hits the next id
                        currLine = celledProjFile{i};
                        
                        while i<length(celledProjFile) && length(currLine)>4 && (linksFlag || ~strncmpi(currLine,'<ID>',4))
                            if (strncmpi(currLine,'<InLinks>',length('<InLinks>')) ...
                                    || strncmpi(currLine,'<OutLinks>',length('<OutLinks>')) ...
                                    || strncmpi(currLine,'<CaActivation>',length('<CaActivation>'))...
                                    || strncmpi(currLine,'<CaDeactivation>',length('<CaDeactivation>'))...
                                    || strncmpi(currLine,'<Gain>',length('<Gain>'))...
                                    || strncmpi(currLine,'<LengthTension>',length('<LengthTension>'))...
                                    || strncmpi(currLine,'<StimulusTension>',length('<StimulusTension>')))
                                linksFlag = 1;
                            end
                            
                            if (strncmpi(currLine,'<CaActivation>',length('<CaActivation>')))
                                caA_m_flag = 1;
                            elseif (strncmpi(currLine,'<CaDeactivation>',length('<CaDeactivation>')))
                                caD_h_flag = 1;
                            elseif (strncmpi(currLine,'</CaActivation>',length('</CaActivation>')))
                                caA_m_flag = 0;
                            elseif (strncmpi(currLine,'</CaDeactivation>',length('</CaDeactivation>')))
                                caD_h_flag = 0;
                            end
                            if (strncmpi(currLine,'</InLinks>',length('</InLinks>')) ...
                                    || strncmpi(currLine,'</OutLinks>',length('</OutLinks>'))...
                                    || strncmpi(currLine,'</CaActivation>',length('</CaActivation>'))...
                                    || strncmpi(currLine,'</CaDeactivation>',length('</CaDeactivation>'))...
                                    || strncmpi(currLine,'</Gain>',length('</Gain>'))...
                                    || strncmpi(currLine,'</LengthTension>',length('</LengthTension>'))...
                                    || strncmpi(currLine,'</StimulusTension>',length('</StimulusTension>')))
                                linksFlag = 0;
                            end
                            
                            if length(currLine)>length(propertyName)+length('</>') && strncmpi(currLine,['<' propertyName],length(propertyName)+1)
                                if skipEntry
                                    %do nothing
                                    
                                    disp('Doing nothing.')
                                    skipEntry = 0;
                                else
                                    
                                    % if the property is found, read/write the property value
                                    compLine3 = currLine(length(currLine)-1:length(currLine));
                                    % swap out the info in the cell array
                                    if ((mModFlag && caA_m_flag) || (hModFlag && caD_h_flag) || ...
                                            (~mModFlag && ~hModFlag && ~caA_m_flag && ~caD_h_flag))
                                        if strcmp(compLine3,'/>')
                                            % then the files are embedded in the tag, and look
                                            % for the beginning of the value data field
                                            startIndex = strfind(currLine,'Value="');
                                            startIndex = startIndex(1)+length('Value="');
                                            % (current line up through value=") + value +
                                            % ' Scale="' + scale + ' Actual="' + actual + '"/>'
                                            
                                            newString = [currLine(1:startIndex-1) num2str(value) '" Scale="' scale '" Actual="' num2str(actual) '"/>'];
                                            
                                            celledProjFile{i} = newString;
                                            fileWithEdit = celledProjFile;
                                            replaced = 1;
                                            break;
                                        else
                                            %We have a nonstandard property. These
                                            %are usually booleans, which have a
                                            %different format in the .aproj files
                                            newString = ['<',propertyName,'>',value,'</',propertyName,'>'];
                                            celledProjFile{i} = newString;
                                            fileWithEdit = celledProjFile;
                                            replaced = 1;
                                            break;
                                            
                                        end
                                    end
                                end
                            else
                                %Do nothing
                            end
                            
                            i = i+1;
                            currLine = celledProjFile{i};
                        end
                    end
                end
                i = i + 1;
            end
            if isempty(fileWithEdit)
                warning(['The property "',propertyName,'" may not exist within the simulation. Please add it to the simulation, save it, and export another standalone simulation.']);
            end
        end
    end
end


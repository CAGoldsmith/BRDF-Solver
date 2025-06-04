classdef joint < matlab.mixin.SetGet
    %JOINT An object of a body attached to a proximal body via a hinge
    %joint, actuated by two antagonistic muscles and controlled by the
    %muscle feedback control.
    %   Geometry and dynamical properties are used to construct a physical
    %   representation of the joint. Muscle properties, which ultimately
    %   determine the configuration and dynamics of the joint, are
    %   available to be modified, either by the user or an optimization
    %   routine, and will return the steady state position of the limb.
    %   Both position and stiffness are available to be tuned.
    
    properties
        %Static values
        jointNumber %Which joint this is in the leg chain
        legNumber %Which leg this is
        numPoints; %number of positions (also equal to the number of stiffnesses) to train on
        maxAngle;
        minAngle;
        angleRange;
        toPlot; %boolean of whether to plot the residuals or not
        useable; %boolean of whether a joint has muscles and thus is useable or not
        uJoint; %normal vector that describes the axis of the joint, in distal frame
        eulerAnglesDisBody_ProxFrame; %orientation of the distal segment, in proximal frame
        eulerAnglesJoint_ProxFrame; %orientation of the joint, in proximal frame
        possibleRotations; %rotation of the joint
        velocity; %rate of rotation of the joint
        maxTorque; %maximum torque expected on the joint
        C_ProxWRTDis; %rotation of the distal body with respect to the proximal one
        pDisFrame_ProxFrame; %location of the distal frame, in proximal frame
        pJoint_ProxFrame; %location of the joint, in proximal frame
        isCTr; %boolean of whether or not this joint is the CTr joint, and therefore gets tuned differently.
        delEdep = 220;
        positionNetParams;
        jointName; %Name of the joint
        
        Rnet; %R value for the network, that is, acceptable range of membrane voltage.
        
        %Posture and synergy properties
        restingPostureAngle;
        restingPostureActivation;
        stanceDirections; %Foot directions away from "straight backward" at which synergies are tested
        stanceVelocities; %Joint velocities that correspond to the stanceDirections
        bodyTranslation;
        bodyRotation;
        pepAngles; %Joint angles the correspond to body translation and rotation
        aepAngles;
        pepActivations;
        aepActivations;
        G_CpgToInter; %Strength of the synapse from the CPG to the interneuron
        VssCPG; %Steady state voltage of one CPG half-center
        vPositionIn; %Position input neuron's resting voltage, tuned such that the neuron is at rest when theta = 0 (startup).
        iPositionInLowerLimit; %Position input current at joint's lower limit
        iPositionInUpperLimit; %Position input current at joint's upper limit
        swingDuration;
        maxStanceDuration;
        maxDutyCycle;
        translateSynCond;
        rotateSynCond;
        pepRestingPot;
        aepRestingPot;
        extMNtonicStimulus;
        flxMNtonicStimulus;
        gScaleTrans
        gScaleRot

        pExtDis
        pFlxDis
        pJointDis
        pExtProx
        pFlxProx
        
        %Robot specific properties        
        extPositionOutServoC;
        extPositionOutServoD;
        extPositionOutServoMinX;
        extPositionOutServoMaxX;
        extPositionOutServoMinY;
        extPositionOutServoMaxY;
        
        flxPositionOutServoC;
        flxPositionOutServoD;
        flxPositionOutServoMinX;
        flxPositionOutServoMaxX;
        flxPositionOutServoMinY;
        flxPositionOutServoMaxY;
        
        extVelocityOutServoC;
        extVelocityOutServoD;
        extVelocityOutServoMinX;
        extVelocityOutServoMaxX;
        extVelocityOutServoMinY;
        extVelocityOutServoMaxY;
        
        flxVelocityOutServoC;
        flxVelocityOutServoD;
        flxVelocityOutServoMinX;
        flxVelocityOutServoMaxX;
        flxVelocityOutServoMinY;
        flxVelocityOutServoMaxY;
        
        extPositionInServoC;
        extPositionInServoD;
        extPositionInServoMinX;
        extPositionInServoMaxX;
        extPositionInServoMinY;
        extPositionInServoMaxY;
        
        flxPositionInServoC;
        flxPositionInServoD;
        flxPositionInServoMinX;
        flxPositionInServoMaxX;
        flxPositionInServoMinY;
        flxPositionInServoMaxY;
        
        
        extVelocityInServoA; %x offset
        extVelocityInServoB; %amplitude
        extVelocityInServoC; %steepness
        flxVelocityInServoA;
        flxVelocityInServoB;
        flxVelocityInServoC;
        
        fCommAngle; %Function that takes in Uext and Uflx and returns the equilibrium position of the joint
        fCommSpeed; %Function that takes in Uext and Uflx and returns the movement speed of the joint
        
        motorEnabled;
        maxJointVelocity;
        isRobot;
        
        %Optimization and parameter sweeping variables
        parameterMap;
        numParams;
        
        %Simulation parameters
        maxSimulationTime = 7000;
        dt = 1;   
        
    end
    
    methods
        function obj = joint(pExtProx,pExtDis,pFlxProx,pFlxDis,pJointDis,eulerAnglesJoint_ProxFrame,pDisFrame_ProxFrame,eulerAnglesBody,ctr,motorEnabled,springPresent,thetaRange,swingDuration,maxStanceDuration,maxDutyCycle,maxJointVelocity,jointName,isRobot,toPlot,jointNum,legNum,muscleProps,lSpringRest,kPassive,cPassive,nprops,loop_parameters,sprops)
            if nargin >= 21
                obj.toPlot = toPlot;
                obj.jointNumber = jointNum;
                obj.legNumber = legNum;
                obj.pExtDis = pExtDis;
                obj.pFlxDis = pFlxDis;
                obj.pJointDis = pJointDis;
                obj.pExtProx = pExtProx;
                obj.pFlxProx = pFlxProx;
                obj.pDisFrame_ProxFrame = pDisFrame_ProxFrame;
                obj.eulerAnglesDisBody_ProxFrame = eulerAnglesBody;
                obj.eulerAnglesJoint_ProxFrame = eulerAnglesJoint_ProxFrame;
                
                obj.swingDuration = swingDuration;
                obj.maxStanceDuration = maxStanceDuration;
                obj.maxDutyCycle = maxDutyCycle;
                obj.maxJointVelocity = maxJointVelocity;
                if abs(diff(thetaRange)) > 1e-3
                    obj.numPoints = 10;
                    obj.possibleRotations = linspace(thetaRange(1),thetaRange(2),obj.numPoints);
                    obj.maxAngle = thetaRange(2);
                    obj.minAngle = thetaRange(1);
                    obj.angleRange = obj.maxAngle - obj.minAngle;
                else
                    keyboard
                    obj.numPoints = 0;
                    obj.possibleRotations = [];
                    obj.maxAngle = 0;
                    obj.minAngle = 0;
                    obj.angleRange = 0;
                end
                
                if islogical(ctr)
                	obj.isCTr = ctr;
                else
                    obj.isCTr = false;
                    warning('ctr must be a boolean. Setting isCTr to false.')
                end
                
                if islogical(motorEnabled)
                    obj.motorEnabled = motorEnabled;
                else
                    error('Joint input motorEnabled must be a boolean.')
                end
                
                obj.isRobot = isRobot;
                
                %Compute the rotation matrix for the distal body with 0 joint
                %rotation.
                obj.C_ProxWRTDis = obj.Rmat(obj.eulerAnglesDisBody_ProxFrame);
                
                %Compute the rotation axis of the joint. This is given in the
                %distal frame, but should be transformed into the proximal
                %frame.
                obj.uJoint = obj.C_ProxWRTDis*obj.Rmat(obj.eulerAnglesJoint_ProxFrame)*[-1;0;0];
                
                %Calculate the position of the joint in the proximal frame.
                %This is very important because this serves as our new "base
                %point" for rotations.
                obj.pJoint_ProxFrame = obj.pDisFrame_ProxFrame + obj.C_ProxWRTDis*obj.pJointDis;
                
                if obj.isRobot
                    if obj.motorEnabled
                        obj.useable = true;
                    else
                        obj.useable = false;
                    end
                else
                    
                    if isempty(obj.pExtDis) || isempty(obj.pFlxDis) || isempty(obj.pExtProx) || isempty(obj.pFlxProx)
                        obj.useable = 0;
                        %If a joint doesn't have any muscles acting on it, we
                        %cannot apply a torque and thus it isn't useable.
                        obj.pExtDis = NaN(3,1);
                        obj.pFlxDis = NaN(3,1);
                        obj.pExtDis_ProxFrame = NaN(3,1);
                        obj.pFlxDis_ProxFrame = NaN(3,1);
                        obj.pExtProx = NaN(3,1);
                        obj.pFlxProx = NaN(3,1);
                    else
                        obj.useable = 1;
                        %We also need to put the distal joint's attachments into the
                        %proximal frame.
                        obj.pExtDis_ProxFrame = repmat(obj.pDisFrame_ProxFrame,1,size(obj.pExtDis,2)) + obj.C_ProxWRTDis*obj.pExtDis;
                        obj.pFlxDis_ProxFrame = repmat(obj.pDisFrame_ProxFrame,1,size(obj.pFlxDis,2)) + obj.C_ProxWRTDis*obj.pFlxDis;
                    end
                    
                    %Calculate the lever arms of the muscles in the proximal frame.
                    %I'm not sure if this is entirely true for multiple attachment
                    %points.
                    obj.rExt = obj.C_ProxWRTDis*(obj.pExtDis(:,end) - obj.pJointDis);
                    obj.rFlx = obj.C_ProxWRTDis*(obj.pFlxDis(:,end) - obj.pJointDis);
                    
                    %Compute the resting length of the muscles, that is, their
                    %length when initialized.
                    
                    obj.lExtRest = obj.computeMuscleLengthWithAttachments(obj.pExtDis_ProxFrame,obj.pExtProx);
                    obj.lFlxRest = obj.computeMuscleLengthWithAttachments(obj.pFlxDis_ProxFrame,obj.pFlxProx);
                    
                    if springPresent
                        obj.passiveStiffnessEnabled = true;
                    end
                end
                
                obj.jointName = jointName;
            end
            if nargin > 21
                obj.muscleProps = muscleProps;
                obj.kseExt = muscleProps(2,1);
                obj.kpeExt = muscleProps(2,2);
                obj.cExt = muscleProps(2,3);
                obj.lWidthExt = muscleProps(2,4);
                obj.amplitudeExt = muscleProps(2,5);
                obj.steepnessExt = muscleProps(2,6);
                obj.xOffsetExt = muscleProps(2,7);
                obj.yOffsetExt = muscleProps(2,8);
                obj.kseFlx = muscleProps(1,1);
                obj.kpeFlx = muscleProps(1,2);
                obj.cFlx = muscleProps(1,3);
                obj.lWidthFlx = muscleProps(1,4);
                obj.amplitudeFlx = muscleProps(1,5);
                obj.steepnessFlx = muscleProps(1,6);
                obj.xOffsetFlx = muscleProps(1,7);
                obj.yOffsetFlx = muscleProps(1,8);
            end
            if nargin > 22
                obj.lSpringRest = lSpringRest;
                obj.kPassive = kPassive;
                obj.cPassive = cPassive;
            end
            if nargin > 25
                obj.nprops = nprops;
                obj.numNeurons = size(nprops,1);
                obj.numStates = obj.numNeurons + 4;
                obj.jOmega = 2*pi*logspace(-2,2,obj.numFreq);
                obj.inputNode = loop_parameters{1};
                obj.outputNode = loop_parameters{2};
                obj.mnNodes = loop_parameters{3}; %[FlexorNode ExtensorNode]
                obj.feedbackNode = loop_parameters{4};
                obj.feedbackSyn=loop_parameters{5};
                obj.jointOutputNode = loop_parameters{6};
                obj.connectionMap =loop_parameters{7};
                obj.numSynapses = length(obj.connectionMap(:,1));
                obj.sprops = sprops;
            end
            if nargin < 15
                error('Not enough arguments for Joint.m')
            end
            
            obj.Rnet = 20;
        end
        
        function success = computeMuscleForceWorstCase(obj, load)
            if isempty(load)
                torqueLoad = obj.maxTorque;
            else
                torqueLoad = abs(load);
            end
            
            T_ext_ext = zeros(obj.numPoints,1);
            T_flx_ext = zeros(obj.numPoints,1);
            T_ext_flx = zeros(obj.numPoints,1);
            T_flx_flx = zeros(obj.numPoints,1);
            obj.lExt = zeros(obj.numPoints,1);
            obj.lFlx = zeros(obj.numPoints,1);
            
            for i=1:obj.numPoints
                %We will now solve for the muscle forces required to
                %counter the maximum torque with the maximum joint
                %stiffness. This is a smooth function, so we will use our
                %quasi-newton optimizer. We will constrain it such that
                %muscle forces are greater than 0. A large upper bound will
                %be used since none is not an option.
                
                %First, assume the torque is negative. This will reveal the
                %maximum extensor force.
                fNetT = @(x)obj.netJointTorque(x,'ext',obj.possibleRotations(i),-torqueLoad,false,0);
                T = fzero(fNetT,0);
                if T > 0
                    T_ext_ext(i) = T(1);
                    T_flx_ext(i) = 0;
                else
                    fNetT = @(x)obj.netJointTorque(x,'flx',obj.possibleRotations(i),-torqueLoad,false,0);
                    T = fzero(fNetT,0);
                    T_ext_ext(i) = 0;
                    T_flx_ext(i) = T(1);
                end
                   
                %Save the most recent muscle length to build our map of
                %joint angle to extensor length.
                obj.lExt(i) = obj.lExtCurrent;
                
                %Next, assume the torque is positive. This will reveal the
                %maximum flexor force.
                fNetT = @(x)obj.netJointTorque(x,'flx',obj.possibleRotations(i),torqueLoad,false,0);
                T = fzero(fNetT,0);
                if T < 0
                    T_ext_flx(i) = 0;
                    T_flx_flx(i) = T(1);
                else
                    fNetT = @(x)obj.netJointTorque(x,'flx',obj.possibleRotations(i),torqueLoad,false,0);
                    T = fzero(fNetT,0);
                    T_ext_flx(i) = T(1);
                    T_flx_flx(i) = 0;
                end
                
                %Again, save the most recent muscle length to build our map
                %of joint angle to flexor length.
                obj.lFlx(i) = obj.lFlxCurrent;
            end
            
            obj.minLExt = min(obj.lExt);
            obj.maxLExt = max(obj.lExt);
            obj.minLFlx = min(obj.lFlx);
            obj.maxLFlx = max(obj.lFlx);
            
            if obj.toPlot
                figure
                clf
                hold on
                title([obj.jointName,' muscle forces under maximum extending torque'])
                plot(obj.possibleRotations,T_ext_ext,'Linewidth',2)
                plot(obj.possibleRotations,T_flx_ext,'g','Linewidth',2)
                legend('T_{ext}','T_{flx}')
                drawnow
                hold off
                
                figure
                clf
                hold on
                title([obj.jointName,' muscle forces under maximum flexing torque'])
                plot(obj.possibleRotations,T_ext_flx,'Linewidth',2)
                plot(obj.possibleRotations,T_flx_flx,'g','Linewidth',2)
                legend('T_{ext}','T_{flx}')
                drawnow
                hold off
                
                figure
                clf
                hold on
                title([obj.jointName,' muscle lengths for desired position and stiffness'])
                plot(obj.possibleRotations,obj.lExt,'b','Linewidth',2);
                plot(obj.possibleRotations,obj.lFlx,'g','Linewidth',2);
                legend('L_{ext}','L_{flx}')
                drawnow
                hold off
                
            end
            obj.toPlot = 0;
            
            %To determine if the position and force maps are accurate, we
            %can make sure that the angle-length curve is monotonically
            %decreasing (for the extensor), and that its rate of change is
            %sufficiently large. If the slope approaches 0 or becomes
            %larger than 0, our muscle forces will become enormous and our
            %muscle lengths will not have a unique joint rotation.
            dlextdthetaNorm = diff(obj.lExt)/range(obj.lExt);
            dlflxdthetaNorm = diff(obj.lFlx)/range(obj.lFlx);
            if any(dlextdthetaNorm >= 0)
                error(['Joint angle-extensor length mapping is not one-to-one for joint ',obj.jointName]);
            end
            if any(dlflxdthetaNorm <= 0)
                error(['Joint angle-flexor length mapping is not one-to-one for joint ',obj.jointName]);
            end
            
            lRangeNorm = range(obj.lExt)/min(obj.lExt);
            if any(dlextdthetaNorm > -.05) || lRangeNorm < 0.2
                warning(['The ',obj.jointName,' rotation - extensor length mapping is poorly conditioned. Position mapping will not be accurate, and muscle forces will be disproportionally large. Modify muscle attachment points until this is no longer the case.'])
            end
            
            %A factor of 1/0.8=1.25 is applied to all of the tension
            %calculations to account for the L-T relationship, which we
            %limit to 0.8 of the maximum.
            %We will use the biological trend of l_width = 1/3*l_rest
            %unless this will not give us desireable muscle force over the
            %range of motion.
            obj.lWidthExt = obj.lExtRest/3; %m
            obj.lWidthFlx = obj.lFlxRest/3; %m
            
            [~,extActivationMaxLength] = obj.steadyStateTension(1,max(obj.lExt),'ext');
            [~,flxActivationMaxLength] = obj.steadyStateTension(1,max(obj.lFlx),'flx');
            
            [~,extActivationMinLength] = obj.steadyStateTension(1,min(obj.lExt),'ext');
            [~,flxActivationMinLength] = obj.steadyStateTension(1,min(obj.lFlx),'flx');
            
            if extActivationMaxLength < 4/5 || extActivationMinLength < 4/5
                if extActivationMaxLength < extActivationMinLength
                    %ext_activation_max_L is the more constricting case
                    obj.lWidthExt = sqrt(5)*(max(obj.lExt) - obj.lExtRest);
                else
                    obj.lWidthExt = -sqrt(5)*(min(obj.lExt) - obj.lExtRest);
                end
            end
            
            if flxActivationMaxLength < 4/5 || flxActivationMinLength < 4/5
                if flxActivationMaxLength < flxActivationMinLength
                    %flx_activation_max_L is the more constricting case
                    obj.lWidthFlx = sqrt(5)*(max(obj.lFlx) - obj.lFlxRest);
                else
                    obj.lWidthFlx = -sqrt(5)*(min(obj.lFlx) - obj.lFlxRest);
                end
            end
            success = 1;
        end
        
        function netTorque = netJointTorqueSpringDesign(obj,thetaDesired,torqueLoad,kSpring)
            %This function is used to design the passive springs in the
            %joint. Given a torque load, desired theta, and spring constant
            %k, this function will calculate the net torque of the joint.
            %It is developed to vary k_spring and seek output net torque =
            %0. Knowing the theta for standing and the torque req'd to
            %stand. L Spring_rest should already be set up to be at rest
            %when standing. We add a small amount to theta_des so that the
            %spring will produce a force.
            thetaTol = .1;
            if torqueLoad <0
                thetaTol = -thetaTol;
            end
            [l,u] = obj.computeMuscleLength(thetaDesired+thetaTol,'flx');
            rRotated = obj.axisAngleRotation(thetaDesired+thetaTol)*obj.rFlx;
            %Tension is positive. If the muscle is longer than its
            %resting length, T_spring > 0
            T_spring = kSpring*(l - obj.lSpringRest);
            
            musc_torque = obj.uJoint'*cross(rRotated,T_spring*u);
            
            % fprintf('Spring contribution: %f%%.\n',100*obj.uJoint'*cross(r_rotated,T_spring*u)/musc_torque);
            
            
            %Find the net torque
            netTorque = obj.uJoint'*(obj.uJoint*torqueLoad) + musc_torque;
            
        end
        
        function objective = netJointTorque(obj,muscForce,muscle,thetaDesired,torqueLoad,normalizedNetTorque,xDotSpring)
            %This function calculates the residual torque on a joint for
            %the given input forces, desired equilibrium position, and
            %external torque. This can be used for a variety of purposes.
            %For a given angle, the forces required to counter an applied
            %torque can be found by seeking output = 0. For given forces,
            %the equilibrium position can be found again by seeking output
            %= 0. The forces are constrained to provide the maximum
            %stiffness the joint should be able to generate.
            if ischar(muscle)
                muscle = {muscle};
            elseif iscell(muscle)
                %do nothing
            else
                disp('Invalid input for Joint.netJointTorque()')
            end
            
            T = muscForce;
            if any(size(T) ~= [length(muscle),1])
                warning('Invalid force input for Joint.netJointTorque()')
            end
            
            muscTorque = 0;
            
            for i=1:length(muscle)
                
                if strcmp(muscle{i},'ext')
                    %Find the length, l, and direction vector, u, of the
                    %extensor muscle.
                    [l,u] = obj.computeMuscleLength(thetaDesired,'ext');
                    obj.lExtCurrent = l;
                    
                    %Transform the muscle's lever arm to its new position.
                    rRotated = obj.axisAngleRotation(thetaDesired)*obj.rExt;
                    
                elseif strcmp(muscle{i},'flx')
                    [l,u] = obj.computeMuscleLength(thetaDesired,'flx');
                    obj.lFlxCurrent = l;
                    
                    rRotated = obj.axisAngleRotation(thetaDesired)*obj.rFlx;
                end
                
                %Add this muscle's contribution to the joint torque.
                muscTorque = muscTorque + obj.uJoint'*cross(rRotated,T(i)*u);
            end
            
            if obj.passiveStiffnessEnabled
                [l,u] = obj.computeMuscleLength(thetaDesired,'flx');
                rRotated = obj.axisAngleRotation(thetaDesired)*obj.rFlx;
                %Tension is positive. If the muscle is longer than its
                %resting length, T_spring > 0
                T_spring = obj.kPassive*(l - obj.lSpringRest)+obj.cPassive*xDotSpring;
                
                muscTorque = muscTorque + obj.uJoint'*cross(rRotated,T_spring*u);
                % fprintf('Spring contribution: %f%%.\n',100*obj.uJoint'*cross(r_rotated,T_spring*u)/musc_torque);
            end
            
            %Find the net torque
            netTorque = torqueLoad + muscTorque;
            if normalizedNetTorque
                objective = netTorque^2;
            else
                objective = netTorque;
            end
        end
        
        function [length,u] = computeMuscleLength(obj,jointAngle,muscle)
            %Takes in a joint angle and the muscle designation.
            %Returns the length of the muscle, and the unit vector direction
            %of the muscle from attachment to attachment in the local frame
            %(for torque calculation).
            
            if strcmp(muscle,'ext')
                %Compute attachment points.
                try pDisProxCurrent = repmat(obj.pJoint_ProxFrame,1,size(obj.pExtDis,2)) + obj.axisAngleRotation(jointAngle)*(obj.pExtDis_ProxFrame - repmat(obj.pJoint_ProxFrame,1,size(obj.pExtDis_ProxFrame,2)));
                catch
                    error(['Muscle attachments for the ',obj.jointName,' extensor are not properly defined. Reinitialize the object and correctly enter information about muscle attachments.'])
                end
                [length,u] = obj.computeMuscleLengthWithAttachments(pDisProxCurrent,obj.pExtProx);
                obj.lExtCurrent = length;
                
            elseif strcmp(muscle,'flx')
                try pDisProxCurrent = repmat(obj.pJoint_ProxFrame,1,size(obj.pFlxDis,2)) + obj.axisAngleRotation(jointAngle)*(obj.pFlxDis_ProxFrame - repmat(obj.pJoint_ProxFrame,1,size(obj.pFlxDis_ProxFrame,2)));
                catch
                    error(['Muscle attachments for the ',obj.jointName,' flexor are not properly defined. Reinitialize the object and correctly enter information about muscle attachments.'])
                end
                [length,u] = obj.computeMuscleLengthWithAttachments(pDisProxCurrent,obj.pFlxProx);
                obj.lFlxCurrent = length;
                
            elseif strcmp(muscle,'all')
                
                length = zeros(2,1);
                u = zeros(2,1);
                
                pDisProxCurrent = repmat(obj.pJoint_ProxFrame,1,size(obj.pExtDis,2)) + obj.axisAngleRotation(jointAngle)*(obj.pExtDis_ProxFrame - repmat(obj.pJoint_ProxFrame,1,size(obj.pExtDis_ProxFrame,2)));
                [length(1),u(1)] = obj.computeMuscleLengthWithAttachments(pDisProxCurrent,obj.pExtProx);
                obj.lExtCurrent = length;
                
                pDisProxCurrent = repmat(obj.pJoint_ProxFrame,1,size(obj.pFlxDis,2)) + obj.axisAngleRotation(jointAngle)*(obj.pFlxDis_ProxFrame - repmat(obj.pJoint_ProxFrame,1,size(obj.pFlxDis_ProxFrame,2)));
                [length(2),u(2)] = obj.computeMuscleLengthWithAttachments(pDisProxCurrent,obj.pFlxProx);
                obj.lFlxCurrent = length;
                
            end
            
        end
        
        function [tension,actualActivation] = steadyStateTension(obj,act,len,muscle)
            if strcmp(muscle,'ext')
                actualActivation = act*(1-(len - obj.lExtRest)^2/obj.lWidthExt^2);
                tension = obj.kseExt/(obj.kpeExt + obj.kseExt)*(actualActivation + obj.kpeExt*(len-obj.lExtRest));
            elseif strcmp(muscle,'flx')
                actualActivation = act*(1-(len - obj.lFlxRest)^2/obj.lWidthFlx^2);
                tension = obj.kseFlx/(obj.kpeFlx + obj.kseFlx)*(actualActivation + obj.kpeFlx*(len-obj.lFlxRest));
            else
                tension = -1;
            end
        end
        
        function residual = trainPositionController(obj)
            if isempty(obj.lExt) || all(obj.lExt)==0
                %If obj.lExt wasn't previously defined, we need to define
                %it now.
                obj.lExt = zeros(obj.numPoints,1);
                for i = 1:obj.numPoints
                    l = obj.computeMuscleLength(obj.possibleRotations(i),'ext');
                    obj.lExt(i) = l;
                end
            end
            obj.vPositionIn = -60e-3 + 20e-3*((0-min(obj.minAngle))/obj.angleRange);
            output = 1e3*(-.060 - obj.vPositionIn) + linspace(0,20,10)';
            
            A = [obj.lExt.^3,obj.lExt.^2,obj.lExt.^1,obj.lExt.^0];
            b = output;
            
            x = linsolve(A,b);
            
            y = x(1)*obj.lExt.^3+x(2)*obj.lExt.^2+x(3)*obj.lExt+x(4)*obj.lExt.^0;
            
            residual = norm(y-b);
            
            figure
            plot(obj.lExt,output)
            hold on
            plot(obj.lExt,x(1)*obj.lExt.^3+x(2)*obj.lExt.^2+x(3)*obj.lExt+x(4)*obj.lExt.^0,'r+')
            grid on
            hold off
            
            obj.AExtLength = x(1)*1e-9;
            obj.BExtLength = x(2)*1e-9;
            obj.CExtLength = x(3)*1e-9;
            obj.DExtLength = x(4)*1e-9;
            
            obj.iPositionInLowerLimit = max(output)*1e-9;
            obj.iPositionInUpperLimit = min(output)*1e-9;
            
        end
        
        function C = axisAngleRotation(obj, angle)
            c = cos(angle);
            s = sin(angle);
            a1 = obj.uJoint(1);
            a2 = obj.uJoint(2);
            a3 = obj.uJoint(3);
            
            C = [c+a1^2*(1-c), a1*a2*(1-c)-a3*s, a1*a3*(1-c)+a2*s;...
                a1*a2*(1-c)+a3*s, c+a2^2*(1-c), a2*a3*(1-c)-a1*s;...
                a1*a3*(1-c)-a2*s, a2*a3*(1-c)+a1*s, c+a3^2*(1-c)];
            
        end
        
        function Vss = vss(obj,Vr,Gmem,Vpre,Esyn,Gmax,Iext)
            Gsyn = obj.syn(Vpre,-.040,-.060,Gmax);
            Vss = (Gmem*Vr + Gsyn'*Esyn + sum(Iext))/(Gmem + sum(Gsyn));
        end
        
        function obj = setMaxTorque(obj,maxTorque)
            obj.maxTorque = maxTorque;
        end
        
        function success = computeRobotAdapterGains(obj,jointHardwareLimits,saveLoc)
            if isempty(obj.restingPostureAngle)
                error('joint.computeRobotAdapterGains() cannot run because no rest posture has been specified.')
            end
            if ~isempty(jointHardwareLimits{obj.jointNumber,obj.legNumber})
                obj.minAngle = jointHardwareLimits{obj.jointNumber,obj.legNumber}(1);
                obj.maxAngle = jointHardwareLimits{obj.jointNumber,obj.legNumber}(2);
            end
            restingPostureAt = 0;
            
            obj.extMNtonicStimulus = restingPostureAt*obj.Rnet*1e-9; 
            obj.flxMNtonicStimulus = restingPostureAt*obj.Rnet*1e-9; 
            
            m = 1/(1-restingPostureAt);
            b = m-1;
            neutralAngle = m*obj.restingPostureAngle - b;
            
            posROM = obj.maxAngle - neutralAngle;
            negROM = neutralAngle - obj.minAngle;

            maxROM = max(posROM,negROM);
            
            %Find the angular distance from the resting angle to each joint
            %limit.

            n = obj.maxStanceDuration/obj.swingDuration;

            %Set maps from neural activity to servo commanded angle.
            obj.extPositionOutServoC = n*maxROM/(obj.Rnet*1e-3); %rad/volts
            obj.flxPositionOutServoC = -n*maxROM/(obj.Rnet*1e-3);

            obj.extPositionOutServoD = 1/2*obj.restingPostureAngle;
            obj.flxPositionOutServoD = obj.extPositionOutServoD;

            speedFac = 2;
            
            obj.extPositionOutServoMinX = 0;
            obj.flxPositionOutServoMinX = 0;

            obj.extPositionOutServoMaxX = obj.Rnet*1e-3;
            obj.flxPositionOutServoMaxX = obj.Rnet*1e-3;

            obj.extPositionOutServoMinY = obj.extPositionOutServoD;
            obj.flxPositionOutServoMinY = obj.extPositionOutServoD;

            obj.extPositionOutServoMaxY = obj.extPositionOutServoD + obj.extPositionOutServoC*obj.Rnet*1e-3;
            obj.flxPositionOutServoMaxY = obj.flxPositionOutServoD + obj.flxPositionOutServoC*obj.Rnet*1e-3;
            
            %Set maps from servo perceived angle to neural activity
            obj.extPositionInServoC = (obj.Rnet*1e-9)/maxROM;
            obj.flxPositionInServoC = -(obj.Rnet*1e-9)/maxROM;
            
            obj.extPositionInServoD = -obj.extPositionInServoC*obj.restingPostureAngle;
            obj.flxPositionInServoD = -obj.flxPositionInServoC*obj.restingPostureAngle;
            
            obj.extPositionInServoMinX = obj.restingPostureAngle;
            obj.flxPositionInServoMinX = obj.minAngle;
            
            obj.extPositionInServoMaxX = obj.maxAngle;
            obj.flxPositionInServoMaxX = obj.restingPostureAngle;
            
            obj.extPositionInServoMinY = 0;
            obj.flxPositionInServoMinY = obj.Rnet*1e-9;
            
            obj.extPositionInServoMaxY = obj.Rnet*1e-9;
            obj.flxPositionInServoMaxY = 0;
            
            %Set maps from neural activity to servo speed.
            maxSpeed = speedFac*maxROM/obj.swingDuration;
            obj.extVelocityOutServoC = maxSpeed/(obj.Rnet*1e-3);
            obj.flxVelocityOutServoC = maxSpeed/(obj.Rnet*1e-3);

            obj.definefCommSpeed(0.05*obj.Rnet*1e-3,0.05*obj.Rnet*1e-3);
            
            %Set maps from servo perceived speed to neural activity
            obj.extVelocityInServoA = 2*obj.extVelocityOutServoMinY; %8*(12/1024) %x offset %8 servo velocity points, scaled to rad/s
            obj.extVelocityInServoB = (obj.Rnet*1e-9); %amplitude
            obj.extVelocityInServoC = 25*obj.Rnet; %steepness
            
            obj.flxVelocityInServoA = -2*obj.flxVelocityOutServoMinY; %-8*(12/1024); %x offset %8 servo velocity points, scaled to rad/s
            obj.flxVelocityInServoB = (obj.Rnet*1e-9); %amplitude
            obj.flxVelocityInServoC = -25*obj.Rnet; %steepness
            
            h = figure;
            hold on
            sgtitle(strrep(obj.jointName,'_','\_'))
            
            subplot(2,2,1);
            hold on
            u = (-10:.5:30)*1e-3;
            
            extOut = min(max(obj.extPositionOutServoD+obj.extPositionOutServoC*u,obj.extPositionOutServoMinY),obj.extPositionOutServoMaxY);
            
            flxOut = min(max(obj.flxPositionOutServoD+obj.flxPositionOutServoC*u,obj.flxPositionOutServoMaxY),obj.flxPositionOutServoMinY);
            
            obj.fCommAngle = @(Uext,Uflx) min(max(obj.extPositionOutServoD+obj.extPositionOutServoC*Uext,obj.extPositionOutServoMinY),obj.extPositionOutServoMaxY) + min(max(obj.flxPositionOutServoD+obj.flxPositionOutServoC*Uflx,obj.flxPositionOutServoMaxY),obj.flxPositionOutServoMinY);
            
            plot(u,extOut)
            plot(u,flxOut)
            plot(u,obj.fCommAngle(u,zeros(size(u))),':')
            plot(u,obj.fCommAngle(zeros(size(u)),u),':')
            plot(u,neutralAngle+zeros(size(u)),'-.')
            legend('Ext alone','Flx alone','Ext+(Flx=0)','(Ext=0)+Flx','location','east')
            xlabel('u')
            ylabel('\theta_{com}')
            grid on
            
            subplot(2,2,2)
            hold on
            th = linspace(obj.minAngle-obj.angleRange/2,obj.maxAngle+obj.angleRange/2,400);
            extIn = obj.extPositionInServoD+obj.extPositionInServoC*th;
            extIn(th > obj.extPositionInServoMaxX) = obj.extPositionInServoMaxY;
            extIn(th < obj.extPositionInServoMinX) = obj.extPositionInServoMinY;
            flxIn = obj.flxPositionInServoD+obj.flxPositionInServoC*th;
            flxIn(th > obj.flxPositionInServoMaxX) = obj.flxPositionInServoMaxY;
            flxIn(th < obj.flxPositionInServoMinX) = obj.flxPositionInServoMinY;
            plot(th,extIn)
            plot(th,flxIn)
            xlabel('\theta')
            ylabel('u')
            grid on
            
            subplot(2,2,3)
            hold on
            
            
            
            extOut = (u < obj.extVelocityOutServoMinX).*obj.extVelocityOutServoMinY + ...
                (u >= obj.extVelocityOutServoMinX & u <= obj.extVelocityOutServoMaxX).*(obj.extVelocityOutServoD+obj.extVelocityOutServoC*u) + ...
                (u > obj.extVelocityOutServoMaxX).*obj.extVelocityOutServoMaxY;

            
            flxOut = (u < obj.flxVelocityOutServoMinX).*obj.flxVelocityOutServoMinY + ...
                (u >= obj.flxVelocityOutServoMinX & u <= obj.flxVelocityOutServoMaxX).*(obj.flxVelocityOutServoD+obj.flxVelocityOutServoC*u) + ...
                (u > obj.flxVelocityOutServoMaxX).*obj.flxVelocityOutServoMaxY;
                        
            plot(u,extOut)
            plot(u,flxOut)
            plot(u,obj.fCommSpeed(u,zeros(size(u))),':')
            plot(u,obj.fCommSpeed(zeros(size(u)),u),':')
            xlabel('u')
            ylabel('d\theta_{com}/dt')
            grid on
            
            subplot(2,2,4)
            hold on
            thDot = linspace(-maxSpeed,maxSpeed,4000);
            extIn = obj.extVelocityInServoB./(1 + exp(obj.extVelocityInServoC*(obj.extVelocityInServoA-thDot)));
            flxIn = obj.flxVelocityInServoB./(1 + exp(obj.flxVelocityInServoC*(obj.flxVelocityInServoA-thDot)));
            plot(thDot,extIn)
            plot(thDot,flxIn)
            xlabel('d\theta/dt')
            ylabel('u')
            grid on    

            savefig(h,[saveLoc '\' obj.jointName])
            
            success = true;
        end
        
        function success = definefCommSpeed(obj,extVminX,flxVminX)
            obj.extVelocityOutServoMinX = extVminX;
            obj.flxVelocityOutServoMinX = flxVminX;
            
            obj.extVelocityOutServoMinY = obj.extVelocityOutServoC*obj.extVelocityOutServoMinX;
            obj.flxVelocityOutServoMinY = obj.flxVelocityOutServoC*obj.flxVelocityOutServoMinX;
            
            obj.extVelocityOutServoMaxX = 0.75*obj.Rnet*1e-3; %Fast MN functionality
            obj.flxVelocityOutServoMaxX = 0.75*obj.Rnet*1e-3; %Fast MN functionality
%             obj.extVelocityOutServoMaxX = 2*obj.Rnet*1e-3; %No Fast MN functionality
%             obj.flxVelocityOutServoMaxX = 2*obj.Rnet*1e-3; %No Fast MN functionality
            
            obj.extVelocityOutServoD = -obj.extVelocityOutServoMinY;
            obj.flxVelocityOutServoD = -obj.extVelocityOutServoMinY;
            
            obj.extVelocityOutServoMaxY = obj.extVelocityOutServoD + obj.extVelocityOutServoC*obj.extVelocityOutServoMaxX;
            obj.flxVelocityOutServoMaxY = obj.flxVelocityOutServoD + obj.flxVelocityOutServoC*obj.flxVelocityOutServoMaxX;
            
            if obj.extVelocityOutServoMaxY > obj.maxJointVelocity
                obj.extVelocityOutServoMaxY = obj.maxJointVelocity;
                obj.extVelocityOutServoMaxX = (obj.extVelocityOutServoMaxY - obj.extVelocityOutServoD)/obj.extVelocityOutServoC;
            else
                obj.extVelocityOutServoMaxY = obj.maxJointVelocity;
                
            end
            if obj.flxVelocityOutServoMaxY > obj.maxJointVelocity
                obj.flxVelocityOutServoMaxY = obj.maxJointVelocity;
                obj.flxVelocityOutServoMaxX = (obj.flxVelocityOutServoMaxY - obj.flxVelocityOutServoD)/obj.flxVelocityOutServoC;
            else
                obj.flxVelocityOutServoMaxY = obj.maxJointVelocity;
            end
            
            obj.fCommSpeed = @(Uext,Uflx) (Uext < obj.extVelocityOutServoMinX).*obj.extVelocityOutServoMinY + ...
                (Uext >= obj.extVelocityOutServoMinX & Uext <= obj.extVelocityOutServoMaxX).*(obj.extVelocityOutServoD+obj.extVelocityOutServoC*Uext) + ...
                (Uext > obj.extVelocityOutServoMaxX).*obj.extVelocityOutServoMaxY + ...
                (Uflx < obj.flxVelocityOutServoMinX).*obj.flxVelocityOutServoMinY + ...
                (Uflx >= obj.flxVelocityOutServoMinX & Uflx <= obj.flxVelocityOutServoMaxX).*(obj.flxVelocityOutServoD+obj.flxVelocityOutServoC*Uflx) + ...
                (Uflx > obj.flxVelocityOutServoMaxX).*obj.flxVelocityOutServoMaxY;
            
            success = true;
        end
        
        function [reconAngles,bmReachable] = designBodyMotionMap(obj,saveLoc1)
            pepAbsAngles = obj.restingPostureAngle + obj.pepAngles;
            bmReachable = ~isnan(pepAbsAngles);
            
            if abs(min(pepAbsAngles(bmReachable)) - max(pepAbsAngles(bmReachable))) < 1e-6
                noVariation = true;
            else
                noVariation = false;
            end
                        
            fprintf('*** BODY MAPPING SYNAPSE STRENGTHS FOR JOINT %s ***\n',obj.jointName)
            if min(pepAbsAngles(bmReachable)) == obj.minAngle
                disp(['Consider changing the maximum flexion of ',obj.jointName,' from ',num2str(180/pi*obj.possibleRotations(1)),' to ',num2str(180/pi*min(pepAbsAngles(:))),' degrees.'])
            end
            if max(pepAbsAngles(bmReachable)) == obj.maxAngle
                disp(['Consider changing the maximum extension of ',obj.jointName,' from ',num2str(180/pi*obj.possibleRotations(end)),' to ',num2str(180/pi*max(pepAbsAngles(:))),' degrees.'])
            end
            
            
            Elo = 0;
            Erange = obj.Rnet;
            Ehi = Elo + Erange;
            
            %Angles and displacements
            %We design the system to produce the translations and rotations
            %in these vectors. The matrices are the result of meshgridding
            %these vectors.
            transVec = obj.bodyTranslation(1,:);
            transMax = obj.bodyTranslation(1,end);
            transMin = obj.bodyTranslation(1,1);
            rotVec = obj.bodyRotation(:,4);
            rotMax = obj.bodyRotation(end,4);
            rotMin = obj.bodyRotation(1,4);
            
%             vTrans = Elo + Erange*(obj.bodyTranslation - min(transVec))/range(transVec);
%             vRot = Elo + Erange*(obj.bodyRotation - min(rotVec))/range(rotVec);
            vTrans = (Ehi+Elo)/2 + obj.bodyTranslation * Erange/(2*transMax);
            vRot = (Ehi+Elo)/2 + obj.bodyRotation * Erange/(2*rotMax);
            
            %solve the least squares problem to map the body's translation
            %and rotation to the joint's rotation. This is the optimal fit.
            A = [vTrans(bmReachable) - mean(vTrans(bmReachable)),vRot(bmReachable) - mean(vRot(bmReachable))];
            b = obj.pepActivations(bmReachable);
            c = obj.aepActivations(bmReachable);
            
            xPep = linsolve(A,b);
            kPepTrans = xPep(1);
            kPepRot = xPep(2);
            
            xAep = linsolve(A,c);
            kAepTrans = xAep(1);
            kAepRot = xAep(2);
            
            %These relationships should be the same for the AEP and PEP. If
            %they are not, then something weird is happening.
            transErr = abs(abs(kPepTrans) - abs(kAepTrans))/abs(kPepTrans);
            rotErr = abs(abs(kPepRot) - abs(kAepRot))/abs(kPepRot);
            if transErr > 0.01
                warning('%s''s mappings from body translation to AEP and PEP disagree by %3.3f%%. This may cause inaccurate limb positioning.',obj.jointName,100*transErr)
            end
            if rotErr > 0.01
                warning('%s''s mappings from body rotation to AEP and PEP disagree by %3.3f%%. This may cause inaccurate limb positioning.',obj.jointName,100*rotErr)
            end
            
            %What is the synaptic conductivity that corresponds to these x
            %values? Scale our G_copy = 1/17, E_copy = 300.
            %delE = 220
            obj.gScaleTrans = abs(kPepTrans)*obj.Rnet/(obj.delEdep - abs(kPepTrans)*obj.Rnet);
            obj.gScaleRot = abs(kPepRot)*obj.Rnet/(obj.delEdep - abs(kPepRot)*obj.Rnet);
            
            
            
            fprintf('Rotation: %.4g nS\n',obj.gScaleRot*1e3);
            fprintf('Translation: %.4g nS\n',obj.gScaleTrans*1e3);
            
            toInvertTrans = false;
            if kPepTrans > 0
                pepTransDescComm = vTrans;
            else
                disp(['INVERT ',obj.jointName,' TRANSLATION']);
                pepTransDescComm = obj.Rnet-vTrans;
                toInvertTrans = true;
            end
            
            toInvertRot = false;
            if kPepRot > 0
                pepRotDescComm = vRot;
            else %mf(1) < 0 && mf(2) < 0
                disp(['INVERT ',obj.jointName,' ROTATION']);
                pepRotDescComm = obj.Rnet-vRot;
                toInvertRot = true;
            end
            
            pep=abs(obj.pepAngles)*obj.extPositionInServoC*1e9; %MAP TO NEURAL VALUES FIRST            
            
            PEPposBias = -1/2*(obj.gScaleTrans + obj.gScaleRot)*obj.delEdep;
            PEPnegBias = PEPposBias;    
            
            toPlotNet = false;
%             toPlotNet = strcmp(obj.jointName,'LF_FTi');
            [~,mappingErr,obj.translateSynCond,obj.rotateSynCond,obj.pepRestingPot,obj.aepRestingPot,pepNeural,pepPosNeural,pepNegNeural] = ...
                descCommNetwork(kPepTrans,kPepRot,PEPposBias,PEPnegBias,toInvertTrans,toInvertRot,Elo,obj.delEdep,obj.Rnet,vTrans,vRot,pep,toPlotNet);
            
     
            
            reconAngles = pepNeural;
            reconAngles(pepNegNeural > 0) = -reconAngles(pepNegNeural > 0);
            reconAngles = obj.restingPostureAngle + reconAngles/(obj.extPositionInServoC*1e9);
            
            fprintf('PEP+ resting potential: %.4g mV\n',obj.pepRestingPot);
            fprintf('PEP- resting potential: %.4g mV\n\n',obj.aepRestingPot);
                        
%             if noVariation
%                 %There is nothing to change in vPepExt and vPepFlx, because
%                 %they are all 0s.
%                 vPepExt = Erest+zeros(size(vPepExt));
%                 vPepFlx = Erest+zeros(size(vPepExt));
%             else
%                 vPepExt(vPepExt == 0) = min(vPepExt(:)); %-60;
%                 vPepFlx(vPepFlx == 0) = max(max(vPepFlx(vPepFlx ~= 0))); %-60;
%             end
%                         
%             %TURN PEP_EXT AND PEP_FLX INTO ANGLES OF EXTENSION AND FLEXION FOR THE PEP
%             extAngles = obj.CPositionOutServo*vPepExt*1e-3+obj.DPositionOutServo; %?
%             flxAngles = obj.CPositionOutServo*vPepFlx*1e-3+obj.DPositionOutServo; %?
%             
%             reconAngles = pepExt.*extAngles + pepFlx.*flxAngles;
%             
%             noMotion = (pepExt & pepFlx) | (~pepExt & ~pepFlx);
%             reconAngles(noMotion) = (extAngles(noMotion) + flxAngles(noMotion))/2;
%             
%             extErr = abs(extAngles - pepAbsAngles);
%             flxErr = abs(flxAngles - pepAbsAngles);
%             
%             totalErrAngle = min(extErr,flxErr);
%             
%             %Add trans and rot synapses to the properties to write during
%             %apply_all_params
            
            %For each leg, map a map of reachability. Basically, a contour
            %map of the total joint space error as a surface of translation
            %and rotation.
            angleLevels = union(obj.restingPostureAngle:.05:obj.maxAngle,...
            obj.restingPostureAngle:-.05:obj.minAngle);
        
            voltLevels = angleLevels*(1e9*obj.extPositionInServoC) + 1e9*obj.extPositionInServoD;
                        
            f = figure;
            clf
            subplot(2,2,1)
            hold on
            %             surf(obj.bodyTranslation,obj.bodyRotation,map_angles)
            %             colorbar
            contour(obj.bodyTranslation,obj.bodyRotation,pepAbsAngles,angleLevels)
            grid on
            xlabel('translation')
            ylabel('rotation')
            title([obj.jointName,' original data'])
            axis([min(transVec),max(transVec),min(rotVec),max(rotVec)])
            if noVariation
                %do not set the colorbar scale
                caxis([0,1])
            else
            	caxis([min(pepAbsAngles(:)),max(pepAbsAngles(:))])
            end
            hold off
            
            subplot(2,2,2)
            hold on
            %             surf(obj.bodyTranslation,obj.bodyRotation,V_sum);
            %             colorbar
            if noVariation
                %do not set the colorbar scale
                caxis([0,1])
            else
            	caxis([min(pepAbsAngles(:)),max(pepAbsAngles(:))])
            end
            contour(obj.bodyTranslation,obj.bodyRotation,reconAngles,angleLevels)
            %             contour(obj.bodyTranslation,obj.bodyRotation,reconAngles,0,'linewidth',5)
            %             plot3(obj.bodyTranslation(:),obj.bodyRotation(:),pepAbsAngles(:),'b.')
            grid on
            xlabel('translation')
            ylabel('rotation')
            title([obj.jointName,' reconstructed output angles'])
            axis([min(transVec),max(transVec),min(rotVec),max(rotVec)])
            hold off
            
            subplot(2,2,3)
            hold on
            %             surf(obj.bodyTranslation,obj.bodyRotation,recon_angles)
            %             colorbar
            if noVariation
                %do not set the colorbar scale
                caxis([0,1])
            else
                
                caxis([min(pepAbsAngles(:)),max(pepAbsAngles(:))]*(obj.extPositionInServoC*1e9) + (obj.extPositionInServoD*1e9))
            end
            
            contour(obj.bodyTranslation,obj.bodyRotation,pepPosNeural,voltLevels); % WE NEED THIS
            plot3(obj.bodyTranslation(pepPosNeural > Ehi),obj.bodyRotation(pepPosNeural > Ehi),1+pepPosNeural(pepPosNeural > Ehi),'ro')
            plot3(obj.bodyTranslation(pepPosNeural < -Ehi),obj.bodyRotation(pepPosNeural < -Ehi),1+pepPosNeural(pepPosNeural < -Ehi),'go')
            grid on
            xlabel('translation')
            ylabel('rotation')
            title([obj.jointName,' PEP Sum'])
            axis([min(transVec),max(transVec),min(rotVec),max(rotVec)])
            hold off
            
            subplot(2,2,4)
            hold on
            surf(obj.bodyTranslation,obj.bodyRotation,abs(180/pi*mappingErr/(obj.extPositionInServoC*1e9)),'EdgeAlpha',0)
            colorbar
            caxis([0,10])
            
            grid on
            xlabel('translation')
            ylabel('rotation')
            title([obj.jointName,' output error, degrees'])
            axis([min(transVec),max(transVec),min(rotVec),max(rotVec)])
            hold off
            drawnow

            savefig(f,[saveLoc1 '\' obj.jointName])
            close gcf
        end
        
        function success = setCpgParams(obj,G,Vss)
            obj.G_CpgToInter = G;
            obj.VssCPG = Vss;
            success = true;
        end
        
        function success = setRestAngle(obj,angle)
            obj.restingPostureAngle = angle;
            success = true;
        end
        
        function props = returnJointControllerProps(obj)
            
            props{1} = obj.lExt;
            props{2} = obj.lFlx;
            
            props{3} = obj.minLExt;
            props{4} = obj.maxLExt;
            props{5} = obj.minLFlx;
            props{6} = obj.maxLFlx;
            
            props{7} = obj.kseExt;
            props{8} = obj.kseFlx;
            
            props{9} = obj.kpeExt;
            props{10} = obj.kpeFlx;
            
            props{11} = obj.amplitudeExt;
            props{12} = obj.amplitudeFlx;
            
            props{13} = obj.cExt;
            props{14} = obj.cFlx;
            
            props{15} = obj.lWidthExt; %m
            props{16} = obj.lWidthFlx; %m
            
            props{17} = obj.lExtRest;
            props{18} = obj.lFlxRest;
            
            props{19} = obj.steepnessExt;
            props{20} = obj.steepnessFlx;
            
            props{21} = obj.yOffsetExt;
            props{22} = obj.yOffsetFlx;
            
            props{23} = obj.xOffsetExt; %V
            props{24} = obj.xOffsetFlx; %V
            
            props{25} = obj.maxTorque;
            
            props{26} = obj.AExtLength;
            props{27} = obj.BExtLength;
            props{28} = obj.CExtLength;
            props{29} = obj.DExtLength;
            props{30} = obj.CExtVel;
            props{31} = obj.positionNetParams;
        end
        
        function props = returnRobotControllerProps(obj)
            
            props{1} = obj.CPositionInServo;
            props{2} = obj.vPositionIn;
            props{3} = obj.CVelocityOutServo;
            props{4} = obj.DVelocityOutServo;
            props{5} = obj.CPositionOutServo;
            props{6} = obj.DPositionOutServo;
        end
        
        function success = assignJointControllerProps(obj,props)
            
            obj.lExt = props{1};
            obj.lFlx = props{2};
            
            obj.minLExt = props{3};
            obj.maxLExt = props{4};
            obj.minLFlx = props{5};
            obj.maxLFlx = props{6};
            
            obj.kseExt = props{7};
            obj.kseFlx = props{8};
            
            obj.kpeExt = props{9};
            obj.kpeFlx = props{10};
            
            obj.amplitudeExt = props{11};
            obj.amplitudeFlx = props{12};
            
            obj.cExt = props{13};
            obj.cFlx = props{14};
            
            obj.lWidthExt = props{15}; %m
            obj.lWidthFlx = props{16}; %m
            
            obj.lExtRest = props{17};
            obj.lFlxRest = props{18};
            
            obj.steepnessExt = props{19};
            obj.steepnessFlx = props{20};
            
            obj.yOffsetExt = props{21};
            obj.yOffsetFlx = props{22};
            
            obj.xOffsetExt = props{23}; %V
            obj.xOffsetFlx = props{24}; %V
            
            obj.maxTorque = props{25};
            
            obj.AExtLength = props{26};
            obj.BExtLength = props{27};
            obj.CExtLength = props{28};
            obj.DExtLength = props{29};
            obj.CExtVel = props{30};
            
            obj.positionNetParams = props{31};
            
            success = 1;
        end
        
        function jointPropsCell = listAllJointParams(obj,~)
            if isempty(obj.restingPostureActivation) && ~obj.isRobot
                error('Find a static posture before writing parameters.')
            end
            
            if obj.isRobot
                
                propsCell = cell(50,3);
            
                %Object name
                for i=1:7
                    propsCell{i,1} = [obj.jointName,'_Ext_ComPos'];
                end
                for i=8:14
                    propsCell{i,1} = [obj.jointName,'_Flex_ComPos'];
                end
                for i=15:21
                    propsCell{i,1} = [obj.jointName,'_Ext_PercPos'];
                end
                for i=22:28
                    propsCell{i,1} = [obj.jointName,'_Flex_PercPos'];
                end
                
                for i=29:35
                    propsCell{i,1} = [obj.jointName,'_Ext_ComVel'];
                end
                for i=36:42
                    propsCell{i,1} = [obj.jointName,'_Flex_ComVel'];
                end
                for i=43:46
                    propsCell{i,1} = [obj.jointName,'_Ext_PercVel'];
                end
                for i=47:50
                    propsCell{i,1} = [obj.jointName,'_Flex_PercVel'];
                end
                
                if obj.isCTr
                    propsCell{51,1} = [obj.jointName,'_Ext_MN'];
                    propsCell{52,1} = [obj.jointName,'_Flex_MN'];
                    propsCell{53,1} = [obj.jointName,'_PEP+'];
                    propsCell{54,1} = [obj.jointName,'_PEP-'];
                end

                %Object property
                propsCell{1,2} = 'C';
                propsCell{2,2} = 'D';
                propsCell{3,2} = 'UseLimits';
                propsCell{4,2} = 'LowerLimit';
                propsCell{5,2} = 'UpperLimit';
                propsCell{6,2} = 'LowerOutput';
                propsCell{7,2} = 'UpperOutput';
                
                propsCell{8,2} = 'C';
                propsCell{9,2} = 'D';
                propsCell{10,2} = 'UseLimits';
                propsCell{11,2} = 'LowerLimit';
                propsCell{12,2} = 'UpperLimit';
                propsCell{13,2} = 'LowerOutput';
                propsCell{14,2} = 'UpperOutput';
                
                propsCell{15,2} = 'C';
                propsCell{16,2} = 'D';
                propsCell{17,2} = 'UseLimits';
                propsCell{18,2} = 'LowerLimit';
                propsCell{19,2} = 'UpperLimit';
                propsCell{20,2} = 'LowerOutput';
                propsCell{21,2} = 'UpperOutput';
                
                propsCell{22,2} = 'C';
                propsCell{23,2} = 'D';
                propsCell{24,2} = 'UseLimits';
                propsCell{25,2} = 'LowerLimit';
                propsCell{26,2} = 'UpperLimit';
                propsCell{27,2} = 'LowerOutput';
                propsCell{28,2} = 'UpperOutput';
                
                propsCell{29,2} = 'C';
                propsCell{30,2} = 'D';
                propsCell{31,2} = 'UseLimits';
                propsCell{32,2} = 'LowerLimit';
                propsCell{33,2} = 'UpperLimit';
                propsCell{34,2} = 'LowerOutput';
                propsCell{35,2} = 'UpperOutput';
                
                propsCell{36,2} = 'C';
                propsCell{37,2} = 'D';
                propsCell{38,2} = 'UseLimits';
                propsCell{39,2} = 'LowerLimit';
                propsCell{40,2} = 'UpperLimit';
                propsCell{41,2} = 'LowerOutput';
                propsCell{42,2} = 'UpperOutput';
                
                propsCell{43,2} = 'A';
                propsCell{44,2} = 'B';
                propsCell{45,2} = 'C';
                propsCell{46,2} = 'UseLimits';
                
                propsCell{47,2} = 'A';
                propsCell{48,2} = 'B';
                propsCell{49,2} = 'C';
                propsCell{50,2} = 'UseLimits';   
                
                if obj.isCTr
                    propsCell{51,2} = 'TonicStimulus';
                    propsCell{52,2} = 'TonicStimulus';
                    propsCell{53,2} = 'TonicStimulus';
                    propsCell{54,2} = 'TonicStimulus';
                end
                
                %%% POSITION OUTPUTS %%%
                propsCell{1,3} = num2str(obj.extPositionOutServoC);
                propsCell{2,3} = num2str(obj.extPositionOutServoD);
                propsCell{3,3} = 'True';
                propsCell{4,3} = num2str(obj.extPositionOutServoMinX);
                propsCell{5,3} = num2str(obj.extPositionOutServoMaxX);
                propsCell{6,3} = num2str(obj.extPositionOutServoMinY);
                propsCell{7,3} = num2str(obj.extPositionOutServoMaxY);
                
                propsCell{8,3} = num2str(obj.flxPositionOutServoC);
                propsCell{9,3} = num2str(obj.flxPositionOutServoD);
                propsCell{10,3} = 'True';
                propsCell{11,3} = num2str(obj.flxPositionOutServoMinX);
                propsCell{12,3} = num2str(obj.flxPositionOutServoMaxX);
                propsCell{13,3} = num2str(obj.flxPositionOutServoMinY);
                propsCell{14,3} = num2str(obj.flxPositionOutServoMaxY);
                
                %%% POSITION INPUTS %%%
                propsCell{15,3} = num2str(obj.extPositionInServoC);
                propsCell{16,3} = num2str(obj.extPositionInServoD);
                propsCell{17,3} = 'True';
                propsCell{18,3} = num2str(obj.extPositionInServoMinX);
                propsCell{19,3} = num2str(obj.extPositionInServoMaxX);
                propsCell{20,3} = num2str(obj.extPositionInServoMinY);
                propsCell{21,3} = num2str(obj.extPositionInServoMaxY);
                
                propsCell{22,3} = num2str(obj.flxPositionInServoC);
                propsCell{23,3} = num2str(obj.flxPositionInServoD);
                propsCell{24,3} = 'True';
                propsCell{25,3} = num2str(obj.flxPositionInServoMinX);
                propsCell{26,3} = num2str(obj.flxPositionInServoMaxX);
                propsCell{27,3} = num2str(obj.flxPositionInServoMinY);
                propsCell{28,3} = num2str(obj.flxPositionInServoMaxY);
                
                %%% VELOCITY OUTPUTS %%%
                propsCell{29,3} = num2str(obj.extVelocityOutServoC);
                propsCell{30,3} = num2str(obj.extVelocityOutServoD);
                propsCell{31,3} = 'True';
                propsCell{32,3} = num2str(obj.extVelocityOutServoMinX);
                propsCell{33,3} = num2str(obj.extVelocityOutServoMaxX);
                propsCell{34,3} = num2str(obj.extVelocityOutServoMinY);
                propsCell{35,3} = num2str(obj.extVelocityOutServoMaxY);
                
                propsCell{36,3} = num2str(obj.flxVelocityOutServoC);
                propsCell{37,3} = num2str(obj.flxVelocityOutServoD);
                propsCell{38,3} = 'True';
                propsCell{39,3} = num2str(obj.flxVelocityOutServoMinX);
                propsCell{40,3} = num2str(obj.flxVelocityOutServoMaxX);
                propsCell{41,3} = num2str(obj.flxVelocityOutServoMinY);
                propsCell{42,3} = num2str(obj.flxVelocityOutServoMaxY);
                
                %%% VELOCITY INPUTS %%%
                propsCell{43,3} = num2str(obj.extVelocityInServoA);
                propsCell{44,3} = num2str(obj.extVelocityInServoB);
                propsCell{45,3} = num2str(obj.extVelocityInServoC);
                propsCell{46,3} = 'False';
                
                propsCell{47,3} = num2str(obj.flxVelocityInServoA);
                propsCell{48,3} = num2str(obj.flxVelocityInServoB);
                propsCell{49,3} = num2str(obj.flxVelocityInServoC);
                propsCell{50,3} = 'False';
                
                if obj.isCTr
                    propsCell{51,3} = num2str(0);
                    propsCell{52,3} = num2str(0);
                    propsCell{53,3} = num2str(obj.extMNtonicStimulus-obj.pepRestingPot*1e-9);
                    propsCell{54,3} = num2str(0);
                    
                end
                
                if obj.isCTr
                    numProps = 52;
                else
                    numProps = 50;
                end
                
                %Make sure every property has a value
                for i=1:numProps
                    if isempty(propsCell{i,3})
                        warning(['Property ',propsCell{i,1},' ',propsCell{i,2},' has no value.']);
                        propsCell{i,3} = num2str(0);
                    end
                end
            else
                propsCell = cell(31,3);
            
                for i=1:18
                    %extensor_name AND flexor_name PROPERTIES SHOULD BE
                    %ASSIGNED AND USED IN THE FUTURE!!!
                    if i <= 9
                        propsCell{i,1} = [obj.jointName,'-Ext'];
                    else
                        propsCell{i,1} = [obj.jointName,'-Flx'];
                    end
                end

                propsCell{1,2} = 'Kse';
                propsCell{10,2} = 'Kse';
                propsCell{1,3} = num2str(obj.kseExt);
                propsCell{10,3} = num2str(obj.kseFlx);

                propsCell{2,2} = 'Kpe';
                propsCell{11,2} = 'Kpe';
                propsCell{2,3} = num2str(obj.kpeExt);
                propsCell{11,3} = num2str(obj.kpeFlx);

                propsCell{3,2} = 'damping';
                propsCell{12,2} = 'damping';
                propsCell{3,3} = num2str(obj.cExt);
                propsCell{12,3} = num2str(obj.cFlx);

                propsCell{4,2} = 'B';
                propsCell{13,2} = 'B';
                propsCell{4,3} = num2str(obj.amplitudeExt);
                propsCell{13,3} = num2str(obj.amplitudeFlx);

                propsCell{5,2} = 'C';
                propsCell{14,2} = 'C';
                propsCell{5,3} = num2str(obj.steepnessExt*1000);
                propsCell{14,3} = num2str(obj.steepnessFlx*1000);

                propsCell{6,2} = 'A';
                propsCell{15,2} = 'A';
                propsCell{6,3} = num2str(obj.xOffsetExt/1000);
                propsCell{15,3} = num2str(obj.xOffsetFlx/1000);

                propsCell{7,2} = 'D';
                propsCell{16,2} = 'D';
                propsCell{7,3} = num2str(obj.yOffsetExt);
                propsCell{16,3} = num2str(obj.yOffsetFlx);

                propsCell{8,2} = 'Lwidth';
                propsCell{17,2} = 'Lwidth';
                propsCell{8,3} = num2str(obj.lWidthExt);
                propsCell{17,3} = num2str(obj.lWidthFlx);

                propsCell{9,2} = 'RestingLength';
                propsCell{18,2} = 'RestingLength';
                propsCell{9,3} = num2str(obj.lExtRest);
                propsCell{18,3} = num2str(obj.lFlxRest);
                
                propsCell{19,1} = [obj.jointName,' Length adapter'];
                propsCell{19,2} = 'A';
                propsCell{19,3} = num2str(obj.AExtLength);
                
                propsCell{20,1} = [obj.jointName,' Length adapter'];
                propsCell{20,2} = 'B';
                propsCell{20,3} = num2str(obj.BExtLength);

                propsCell{21,1} = [obj.jointName,' Length adapter'];
                propsCell{21,2} = 'C';
                propsCell{21,3} = num2str(obj.CExtLength);

                propsCell{22,1} = [obj.jointName,' Length adapter'];
                propsCell{22,2} = 'D';
                propsCell{22,3} = num2str(obj.DExtLength);
                
                propsCell{23,1} = [obj.jointName,' Length adapter'];
                propsCell{23,2} = 'UseLimits';
                propsCell{23,3} = 'true';
                
                propsCell{24,1} = [obj.jointName,' Length adapter'];
                propsCell{24,2} = 'LowerLimit';
                propsCell{24,3} = num2str(min(obj.lExt));
                
                propsCell{25,1} = [obj.jointName,' Length adapter'];
                propsCell{25,2} = 'UpperLimit';
                propsCell{25,3} = num2str(max(obj.lExt));
                
                propsCell{26,1} = [obj.jointName,' Length adapter'];
                propsCell{26,2} = 'LowerOutput';
                propsCell{26,3} = num2str(obj.iPositionInLowerLimit);
                
                propsCell{27,1} = [obj.jointName,' Length adapter'];
                propsCell{27,2} = 'UpperOutput';
                propsCell{27,3} = num2str(obj.iPositionInUpperLimit);

                propsCell{28,1} = [obj.jointName,' Ext adapter'];
                propsCell{28,2} = 'D';
                propsCell{28,3} = num2str(obj.DMNExt/1000);

                propsCell{29,1} = [obj.jointName,' Flx adapter'];
                propsCell{29,2} = 'D';
                propsCell{29,3} = num2str(obj.DMNFlx/1000);
                
                propsCell{30,1} = [obj.jointName,' Resting Position'];
                propsCell{30,2} = 'RestingPot';
                propsCell{30,3} = num2str(obj.restingPostureActivation/1000);
                
                propsCell{31,1} = [obj.jointName,' Perceived Position'];
                propsCell{31,2} = 'RestingPot';
                propsCell{31,3} = num2str(obj.vPositionIn);
                
                propsCell{32,1} = [obj.jointName, '-Spring'];
                propsCell{32,2} = 'NaturalLength';
                propsCell{32,3} = num2str(obj.lSpringRest);
                
                propsCell{33,1} = [obj.jointName, '-Spring'];
                propsCell{33,2} = 'Stiffness';
                propsCell{33,3} = num2str(obj.kPassive);
                
                propsCell{34,1} = [obj.jointName, '-Spring'];
                propsCell{34,2} = 'Damping';
                propsCell{34,3} = num2str(obj.cPassive);

                for i=1:31
                    if isempty(propsCell{i,3})
                        propsCell{i,3} = num2str(0);
                    end
                end
            end
            
            posturePropsCell = cell(10,3);
            
            %Object name
            posturePropsCell{1,1} = [obj.jointName,'_PEP+'];
            posturePropsCell{2,1} = [obj.jointName,'_PEP-'];
            for i=3:6
                posturePropsCell{i,1} = ['Desc ', obj.jointName,' Trans'];
            end
            for i=7:10
                posturePropsCell{i,1} = ['Desc ', obj.jointName,' Rotate'];
            end
            
            %SYNAPSE STRENGTHS FOR DIRECTED LOCOMOTION CONTROL    
            %Object property
            posturePropsCell{1,2} = 'RestingPot';
            posturePropsCell{2,2} = 'RestingPot';
            posturePropsCell{3,2} = 'Equil';
            posturePropsCell{4,2} = 'SynAmp';
            posturePropsCell{5,2} = 'ThreshV';
            posturePropsCell{6,2} = 'SaturateV';
            posturePropsCell{7,2} = 'Equil';
            posturePropsCell{8,2} = 'SynAmp';
            posturePropsCell{9,2} = 'ThreshV';
            posturePropsCell{10,2} = 'SaturateV';
            
            %Property value
            posturePropsCell{1,3} = num2str(obj.pepRestingPot*1e-3);
            posturePropsCell{2,3} = num2str(obj.aepRestingPot*1e-3);
            
            posturePropsCell{3,3} = num2str(obj.delEdep*1e-3);
            posturePropsCell{4,3} = num2str(obj.translateSynCond*1e-6);
            posturePropsCell{5,3} = '0';
            posturePropsCell{6,3} = '0.020';

            posturePropsCell{7,3} = num2str(obj.delEdep*1e-3);
            posturePropsCell{8,3} = num2str(obj.rotateSynCond*1e-6);
            posturePropsCell{9,3} = '0';
            posturePropsCell{10,3} = '0.020';
                
            jointPropsCell = [propsCell;posturePropsCell];
        end
        
       
        function success = setLocomotionAngles(obj,pepAngles,aepAngles)
            obj.pepAngles = pepAngles';
            try obj.pepActivations = pepAngles'*obj.extPositionInServoC*1e9;
            catch
                keyboard
            end
            obj.aepAngles = aepAngles';
            obj.aepActivations = aepAngles'*obj.extPositionInServoC*1e9;
            success = true;
            
%             obj.definefCommSpeed(0,0);
%             
%             obj.pepActivations = zeros(size(obj.pepAngles));
%             obj.aepActivations = zeros(size(obj.pepAngles));
%             
%             for i=1:length(obj.pepActivations(:))
%                 %PEP Angles
%                 if obj.pepAngles(i) > 0
%                     if obj.isCTr
% %                     f = @(x) obj.fCommAngle(x,0) - obj.restingPostureAngle - obj.pepAngles(i);
%                         keyboard
%                     else
%                         f = @(x) obj.fCommSpeed(x,0) - 2*obj.pepAngles(i)/obj.swingDuration;
%                         lims = [obj.extVelocityOutServoMinX,obj.extVelocityOutServoMaxX];
%                         try obj.pepActivations(i) = fzero(f,lims);
%                         catch
%                             obj.pepActivations(i) = NaN;
%                         end
%                     end
%                 elseif obj.pepAngles(i) < 0
%                     if obj.isCTr
% %                     f = @(x) obj.fCommAngle(0,x) - obj.restingPostureAngle - obj.pepAngles(i);
%                         keyboard
%                     else
%                         f = @(x) obj.fCommSpeed(0,x) + 2*obj.pepAngles(i)/obj.swingDuration;
%                         lims = [obj.flxVelocityOutServoMinX,obj.flxVelocityOutServoMaxX];
%                         try obj.pepActivations(i) = -fzero(f,lims);
%                         catch
%                             obj.pepActivations(i) = NaN;
%                         end
%                     end
%                 else
%                     obj.pepActivations(i) = 0;
%                 end
%                 %AEP Angles
%                 if obj.aepAngles(i) > 0
%                     if obj.isCTr
% %                     f = @(x) obj.fCommAngle(x,0) - obj.restingPostureAngle - obj.aepAngles(i);
%                         keyboard
%                     else
%                         f = @(x) obj.fCommSpeed(x,0) - 2*obj.aepAngles(i)/obj.swingDuration;
%                         lims = [obj.extVelocityOutServoMinX,obj.extVelocityOutServoMaxX];
%                         try obj.aepActivations(i) = fzero(f,lims);
%                         catch
%                             obj.aepActivations(i) = NaN;
%                         end
%                     end
%                 elseif obj.aepAngles(i) < 0
%                     if obj.isCTr
% %                     f = @(x) obj.fCommAngle(0,x) - obj.restingPostureAngle - obj.aepAngles(i);
%                         keyboard
%                     else
%                         f = @(x) obj.fCommSpeed(0,x) + 2*obj.aepAngles(i)/obj.swingDuration;
%                         lims = [obj.flxVelocityOutServoMinX,obj.flxVelocityOutServoMaxX];
%                         try obj.aepActivations(i) = -fzero(f,lims);
%                         catch
%                             obj.aepActivations(i) = NaN;
%                         end
%                     end
%                 else
%                     obj.aepActivations(i) = 0;
%                 end
%                 
%             end
            
            
            tempPepActivations = obj.pepActivations;
            
            th = 1;
            
            f = @(thresh) 1 - numel(tempPepActivations(abs(tempPepActivations) < thresh))./numel(tempPepActivations);
            while f(th) < 0.9 && th > 0
                th = th - 2e-2;
            end
            
            obj.definefCommSpeed(th*1e-3,th*1e-3);
            
            fprintf('For joint %s, th became %0.3f and the velocity lower limit is %0.3f.\n',obj.jointName,th,obj.extVelocityOutServoMinX);
            
%             obj.pepActivations = 1e3*obj.pepActivations;
%             obj.aepActivations = 1e3*obj.aepActivations;
%                    
%             if obj.isCTr
%                 keyboard
%             end
            
        end        
    end
    
    methods(Static)
        function [network_params,solution_error] = trainHiddenLayer(input,output,n_hid_nodes_vec,error_thresh,GA_pop_size,to_plot_train,figure_title)
            %How many network configurations (numbers of hidden nodes) will
            %we try?
            num_trials = length(n_hid_nodes_vec);
            
            %Sort the nodes vector, to make sure we are trying increasing
            %numbers of nodes. This increases the chances that we'll find
            %the simplest possible network for the desired mapping, and
            %will also hopefully shorten training time (since we try simple
            %networks first, and will stop if an acceptable solution is
            %found).
            n_hid_nodes_vec = sort(n_hid_nodes_vec);
            
            %The user may specify a different GA population size for each
            %configuration, or a single value to be used for all. 5000 is a
            %good value.
            if length(GA_pop_size) == 1 && num_trials > 1
                GA_pop_size = zeros(num_trials,1) + GA_pop_size;
            end
            
            %How many boundary nodes (input + output) do we have? This will
            %set the number of synapses.
            n_bound_nodes = size(input,2) + size(output,2);
            
            solutions = NaN(3*n_bound_nodes*max(n_hid_nodes_vec),num_trials);
            refined_solutions = solutions;
            
            errors = NaN(num_trials,1);
            
            add_another_node = 1;
            i = 1;
            
            while i <= num_trials && add_another_node
                n_hid_nodes = int32(n_hid_nodes_vec(i));
                n_syn = n_bound_nodes*n_hid_nodes;
                f = @(x)network_fitness_mex(x,input,output,n_hid_nodes);
                con_E = repmat([-100,0],n_syn,1); %Bound on the synapse potential
                con_L = repmat([-80,-40],n_syn,1); %Bound on the lower voltage
                con_R = repmat([1,20],n_syn,1); %Bound on the voltage range
                bnds = [con_E;con_L;con_R];
                
                pop_size = GA_pop_size(i);
                
                %Compute GA solution, and network output
                max_gen = 50;
                percent_same_for_convergence = 0.95;
                radius_for_convergence = 1e-6;
                mating_fraction = 0.5;
                mutation_prob = .001;
                random_seed = 1;
                to_parallelize = false;
                to_print = 1;
                
                params = GA(f,bnds,[],[max_gen,percent_same_for_convergence,...
                    radius_for_convergence],pop_size,mating_fraction,...
                    'single',mutation_prob,random_seed,to_parallelize,...
                    to_plot_train,to_print);
                
                solutions(1:3*n_syn,i) = params;
                [~,~] = network_fitness(params,input,output,n_hid_nodes);
                
                %Refine the GA solution, and network output
                %                 refined_params = bnd_nelder_mead(f,params,bnds,[],1e4,[],[],to_print);
                refined_params = bnd_con_pen_nelder_mead(f,params,bnds,[],[],1e4,[],[],to_print);
                [error,~] = network_fitness(refined_params,input,output,n_hid_nodes);
                refined_solutions(1:3*n_syn,i) = refined_params;
                
                if error < error_thresh
                    add_another_node = 0;
                end
                
                errors(i) = error;
                i = i + 1;
            end
            
            [~,ind] = min(errors);
            
            [~,GA_performance] = network_fitness(solutions(1:3*n_bound_nodes*n_hid_nodes_vec(ind),ind),input,output,ind);
            [~,refined_performance] = network_fitness(refined_solutions(1:3*n_bound_nodes*n_hid_nodes_vec(ind),ind),input,output,ind);
            
            to_plot_train = 1;
            if to_plot_train
                
                if n_bound_nodes == 2
                    figure
                    clf
                    hold on
                    titlestr = strcat(figure_title,' (',num2str(ind),' hidden nodes)');
                    title(titlestr)
                    plot(input,output,'+','Linewidth',2)
                    plot(input,GA_performance(:,1),'g.','Linewidth',2)
                    plot(input,refined_performance(:,1),'ro','Linewidth',2)
                    legend('Desired','GA Result','NM Refined')
                    xlabel('Input Voltage (mV)')
                    ylabel('Output Voltage (mV)')
                    drawnow
                    hold off
                elseif n_bound_nodes == 3
                    figure
                    clf
                    hold on
                    titlestr = strcat(figure_title,' (',num2str(ind),' hidden nodes)');
                    title(titlestr)
                    plot3(input(:,1),input(:,2),output,'+','Linewidth',2)
                    plot3(input(:,1),input(:,2),GA_performance(:,1),'g.','Linewidth',2)
                    plot3(input(:,1),input(:,2),refined_performance(:,1),'ro','Linewidth',2)
                    legend('Desired','GA Result','NM Refined')
                    xlabel('Input 1 Voltage (mV)')
                    ylabel('Input 2 Voltage (mV)')
                    zlabel('Output Voltage (mV)')
                    drawnow
                    hold off
                else
                    %Do nothing
                    disp('Network mapping is too high dimensionality to plot.')
                end
            end
            
            %remove the NaN characters from the vector of final values.
            %E,V_lo,V_hi
            temp_ref_sol = reshape(refined_solutions(1:3*n_bound_nodes*ind,ind),[],3);
            temp_ref_sol(:,3) = temp_ref_sol(:,2) + temp_ref_sol(:,3);
            network_params = temp_ref_sol;
            solution_error = errors(ind);
            %             fprintf('Best number of nodes: %i\n',ind);
        end
        
        function G = syn(V,hi,lo,Gmax)
            G = Gmax.*min(max((V-lo)/(hi-lo),0),1);
        end
        
        function [musc_len,f_vec] = computeMuscleLengthWithAttachments(dist,prox)
            %Inputs are the distal and proximal muscle attachments.
            %Find the length of the muscle and the direction of its force.
            %This function is intended to be called by other functions
            %within this class; it is not for the "user."
            
            %First we need to reverse the order of the proximal
            %attachments. This way we are calculating the distance between
            %attachments in the proper order.
            order = size(prox,2):-1:1;
            attachment_pts = [dist,prox(:,order)];
            
            musc_len = 0;
            
            %Add the magnitude of the vectors between each attachment
            for i=1:size(attachment_pts,2)-1
                musc_len = musc_len + norm(attachment_pts(:,i) - attachment_pts(:,i+1));
            end
            
            %Since in Leg we ask the user to list the muscle attachments
            %from the terminus toward the other segment, we should use the
            %last proximal and last distal attachments to determine a unit
            %vector along which the muscle force acts. It should point
            %toward the proximal body, so the proximal point is the head of
            %the vector
            f_vec = (prox(:,end) - dist(:,end))/norm(prox(:,end) - dist(:,end));
            
        end
        
        function Cmat = Rmat(angs)
            %Take in a euler angle triplet, and return the associated
            %rotation matrix.
            c1 = cos(angs(1));
            s1 = sin(angs(1));
            c2 = cos(angs(2));
            s2 = sin(angs(2));
            c3 = cos(angs(3));
            s3 = sin(angs(3));
            %XYZ
            Cmat = [c2*c3,-c2*s3,s2;c1*s3+c3*s1*s2,c1*c3-s1*s2*s3,-c2*s1;s1*s3-c1*c3*s2,c3*s1+c1*s2*s3,c1*c2];
        end
        
        
    end
end
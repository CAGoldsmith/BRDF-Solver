classdef leg < matlab.mixin.SetGet
    %JACOBIAN Create a Jacobian manipulator object from an AnimatLab
    %simulation file.
    %   Detailed explanation goes here
    properties
        originalText; %content of the .asim file
        jointTypes; %type of joints in the leg, as read from the .asim file
        numBodies; %number of bodies in the leg, as given by the user's list
        organismName; %name of the organism in the simulation file that we are interested in
        legNumber; %Which leg in the organism is this?
        bodies; %names of the segments in the leg
        mass; %mass of the segments in the leg
        lengths; %length of the segments in the leg
        width; %width of the segments in the leg
        height; %height of the segments in the leg
        plasticL; %length of the plastic portions of the segments of the leg
        plasticW; %width of the plastic portions of the segments of the leg
        plasticH; %height of the plastic portions of the segments of the leg
        plasticMass; %mass of the plastic portions of the segments of the leg
        plasticCOM; %COM of the plastic portions, with reference to the most proximal joint it attaches to
        joints; %names of the joints in the leg
        toPlot; %boolean of whether or not to plot graphs of the leg.
        vecLen; %length of the axis vector drawn in the kinematic maps.
        gst0Bodies; %Initial exponential twists of the bodies
        gst0Joints; %Initial exponential twists of the joints
        gst0COMs; %Initial exponential twists of the centers of mass
        twist; %Joint twists for this leg
        posBodies; %position of the distal body in the proximal frame
        posJoints; %location of the distal joint in the proximal frame
        posAttachments; %location of muscle attachments in their own frame
        jointObj; %cell array of joint objects (for each joint in the leg)
        eulerAnglesJoints;
        projectFile; %project file
        toePositionKnown = 0; %Boolean about whether the toe position is known or not.
        footVec; %vector that describes the toe wrt the most distal joint. Important for calculating forces from the foot.
        ankleFactor; %multiply distal max torque by 3 if it operates like an ankle.
        jointModelLowerLimits; %lower limit of each joint
        jointModelUpperLimits; %upper limit of each joint
        motorEnabled; %keeps track of which joints have motors enabled
        jointHardwareLimits; %limits of motion for each joint
        isRobot; %if isRobot, then muscles are ignored
        aepConfig; %Configuration of the resting pose. This is fully flexed.
        restingConfig; %Configuration when standing.
        lengthScale; %Factor to scale distances to meters.
        swingDuration; %Length of swing phase in seconds.
        maxStanceDuration; %Maximum length of stance phase in seconds.
        maxDutyCycle; %Maximum duty cycle during walking (at lowest speed)
        ctrIndex; %Index of the joint that is the CTr, which, as the lifter, is tuned differently.
        connectionMap;
    end

    methods
        function obj = leg(projectFile,organismName,legNumber,bodies,joints,ctr,swingDuration,maxStanceDuration,maxDutyCycle,maxJointVelocity,isRobot,toPlot,legPlastic,loopParameters)
            %Initialize some variables.
            obj.organismName = organismName;
            obj.jointTypes = cell(1,length(joints));
            obj.isRobot = isRobot;
            obj.ctrIndex = ctr;
            obj.legNumber = legNumber;
            obj.swingDuration = swingDuration;
            obj.maxStanceDuration = maxStanceDuration;
            obj.maxDutyCycle = maxDutyCycle;

            if nargin == 15
                obj.connectionMap = loopParameters{7};
            end

            obj.numBodies = 0;
            i = 1;
            while i <= length(bodies) && ~isempty(bodies{i})
                obj.numBodies = obj.numBodies + 1;
                i = i + 1;
            end

            %Populate the appropriate portions of the plastic properties
            if ~isempty(legPlastic)
                obj.plasticL = legPlastic.length{legNumber};
                obj.plasticW = legPlastic.width{legNumber};
                obj.plasticH = legPlastic.height{legNumber};
                obj.plasticCOM = legPlastic.COM{legNumber};
                obj.plasticMass = legPlastic.mass{legNumber};
            end
            %Load the text of the simulation file
            try obj.originalText = importdata(projectFile);
            catch
                error('No simulation file exists, or the directory is incorrect. Double check the simulation directory, and open AnimatLab and export a standalone simulation.')
            end

            %Find the scaling of the simulation.
            lengthUnitsFound = strfind(obj.originalText,'<DistanceUnits>');
            lengthUnitsInd = find(~cellfun(@isempty,lengthUnitsFound),1,'first');
            lengthUnitsStr = obj.originalText{lengthUnitsInd};

            startInd = strfind(lengthUnitsStr,'<DistanceUnits>')+15;
            endInd = strfind(lengthUnitsStr,'</DistanceUnits>')-1;
            lengthString = lengthUnitsStr(startInd:endInd);

            if strcmp(lengthString,'Millimeters')
                obj.lengthScale = 1e-3;
            elseif strcmp(lengthString,'Centimeters')
                obj.lengthScale = 1e-2;
            elseif strcmp(lengthString,'Decimeters')
                obj.lengthScale = 1e-1;
            elseif strcmp(lengthString,'Meters')
                obj.lengthScale = 1;
            elseif strcmp(lengthString,'Decameters')
                obj.lengthScale = 1e1;
            elseif strcmp(lengthString,'Centameters')
                obj.lengthScale = 1e2;
            elseif strcmp(lengthString,'Kilometers')
                obj.lengthScale = 1e3;
            else
                warning('Length scale is not recognized. Defaulting to meters.')
                obj.lengthScale = 1;
            end

            %Find where the list of organisms begins
            %             organismsFound = strfind(obj.originalText,'<Organisms>');
            %             organismsLowerLimit = find(~cellfun(@isempty,organismsFound),1,'first');
            organismsLowerLimit = find(contains(obj.originalText,'<Organisms>'),1,'first');

            %Find where the organism we want begins
            %             organismFound = strfind(obj.originalText(organismsLowerLimit:end),['<Name>',obj.organismName,'</Name>']);
            %             organismLowerLimit = find(~cellfun(@isempty,organismFound),1,'first') + organismsLowerLimit - 1;
            organismLowerLimit = find(contains(obj.originalText(organismsLowerLimit:end),['<Name>',obj.organismName,'</Name>']),1,'first') + organismsLowerLimit - 1;
            %Rotation of each body in the inertial frame
            obj.gst0Bodies = NaN(4,4,obj.numBodies);
            obj.gst0Joints = NaN(4,4,obj.numBodies);

            %Position of the joints and bodies, in their local frames.
            obj.posJoints = NaN(3,obj.numBodies);
            obj.posBodies = NaN(3,obj.numBodies);
            obj.posAttachments = cell(obj.numBodies,1); %one for each joint
            eulerAnglesBodies = NaN(3,obj.numBodies);
            obj.eulerAnglesJoints = NaN(3,obj.numBodies);

            %Initialize other object values
            obj.bodies = bodies;
            obj.joints = joints;
            obj.toPlot = toPlot;
            springPresent = false(obj.numBodies,1);
            lSpringRest = NaN(1,length(springPresent));
            kPassive = NaN(1,length(springPresent));
            cPassive = NaN(1,length(springPresent));

            %Initialize joint range
            obj.jointModelLowerLimits = NaN(obj.numBodies,1);
            obj.jointModelUpperLimits = NaN(obj.numBodies,1);

            obj.motorEnabled = false(obj.numBodies,1);

            i = 1;
            while i <= obj.numBodies
                %Find the body of interest

                if ~isempty(obj.bodies{i})
                    chainLowerLimit = find(contains(obj.originalText,['<Name>',obj.bodies{i},'</Name>']),1,'first');

                    if isempty(chainLowerLimit)
                        error(['<Name>',obj.bodies{i},'</Name> not found. An improper body or joint was entered, or perhaps out of order.'])
                    end

                    %Find the string that lists its rotation
                    bodyRotationInd = find(contains(obj.originalText(chainLowerLimit:end),'<Rotation'),1,'first') + chainLowerLimit - 1;
                    bodyRotationStr = obj.originalText(bodyRotationInd+1:bodyRotationInd+3);

                    for j=1:3
                        startInd = strfind(bodyRotationStr{j},'Actual="')+8;
                        endInd = strfind(bodyRotationStr{j},'"/>')-1;
                        eulerAnglesBodies(j,i) = pi/180*str2double(bodyRotationStr{j}(startInd:endInd));
                    end

                    %                     %Compute the rotation matrix for each axis.
                    %                     rotx_b = [1,0,0;0,cos(eulerAnglesBodies(1,i)),-sin(eulerAnglesBodies(1,i));0,sin(eulerAnglesBodies(1,i)),cos(eulerAnglesBodies(1,i))];
                    %                     roty_b = [cos(eulerAnglesBodies(2,i)),0,sin(eulerAnglesBodies(2,i));0,1,0;-sin(eulerAnglesBodies(2,i)),0,cos(eulerAnglesBodies(2,i))];
                    %                     rotz_b = [cos(eulerAnglesBodies(3,i)),-sin(eulerAnglesBodies(3,i)),0;sin(eulerAnglesBodies(3,i)),cos(eulerAnglesBodies(3,i)),0;0,0,1];
                    %
                    %Find the position of the body in the local frame, that
                    %is, with respect to the parent body.
                    bodyPositionInd = find(contains(obj.originalText(chainLowerLimit:end),'<LocalPosition'),1,'first') + chainLowerLimit - 1;
                    bodyPositionStr = obj.originalText(bodyPositionInd+1:bodyPositionInd+3);

                    for j=1:3
                        startInd = strfind(bodyPositionStr{j},'Actual="')+8;
                        endInd = strfind(bodyPositionStr{j},'"/>')-1;
                        obj.posBodies(j,i) = str2double(bodyPositionStr{j}(startInd:endInd));
                    end

                    %Find the exponential of the twists of the body
                    %segments
                    %Find the position of the body in the local frame, that
                    %is, with respect to the parent body.
                    bodyTwistInd = find(contains(obj.originalText(chainLowerLimit:end),'<LocalMatrix'),1,'first') + chainLowerLimit - 1;
                    bodyTwistStr = obj.originalText{bodyTwistInd};

                    g = NaN(4);
                    endInds = strfind(bodyTwistStr,',');
                    endInds = [endInds,strfind(bodyTwistStr,'</LocalMatrix>')]; %#ok<AGROW>
                    startInd = strfind(bodyTwistStr,'<LocalMatrix>')+13;
                    for j=1:16
                        endInd = endInds(j)-1;
                        g(j) = str2double(bodyTwistStr(startInd:endInd));
                        startInd = endInd+1;
                    end
                    g(1:3,4) = g(1:3,4)*obj.lengthScale; %convert to meters.
                    obj.gst0Bodies(:,:,i) = g;

                    %keyboard

                    %Find the CoM in the local frame
                    bodyCOMInd = find(contains(obj.originalText(chainLowerLimit:end),'<COM'),1,'first') + chainLowerLimit - 1;

                    bodyCOMxInd = find(contains(obj.originalText(bodyCOMInd:end),'<X Value'),1,'first') + bodyCOMInd - 1;
                    bodyCOMxStr = obj.originalText{bodyCOMxInd}(10:end);

                    startInd = strfind(bodyCOMxStr,'Actual="')+8;
                    endInd = strfind(bodyCOMxStr,'"/>')-1;
                    COMlf(1,1) = str2double(bodyCOMxStr(startInd:endInd));

                    bodyCOMyInd = find(contains(obj.originalText(bodyCOMInd:end),'<Y Value'),1,'first') + bodyCOMInd - 1;
                    bodyCOMyStr = obj.originalText{bodyCOMyInd}(10:end);

                    startInd = strfind(bodyCOMyStr,'Actual="')+8;
                    endInd = strfind(bodyCOMyStr,'"/>')-1;
                    COMlf(2,1) = str2double(bodyCOMyStr(startInd:endInd));

                    bodyCOMzInd = find(contains(obj.originalText(bodyCOMInd:end),'<Z Value'),1,'first') + bodyCOMInd - 1;
                    bodyCOMzStr = obj.originalText{bodyCOMzInd}(10:end);

                    startInd = strfind(bodyCOMzStr,'Actual="')+8;
                    endInd = strfind(bodyCOMzStr,'"/>')-1;
                    COMlf(3,1) = str2double(bodyCOMzStr(startInd:endInd));

                    R = g(1:3,1:3);

                    COMgf = R*COMlf + g(1:3,4);

                    COMg = g(1:4,1:3);
                    COMg(:,4) = [COMgf; 1];


                    obj.gst0COMs(:,:,i) = COMg;

                    %What type of shape is it?
                    bodyShapeInd = find(contains(obj.originalText(chainLowerLimit:end),'<Type>'),1,'first') + chainLowerLimit - 1;
                    bodyShape = obj.originalText{bodyShapeInd}(7:end-7);

                    %Extract the mass for the segment

                    massInd = find(contains(obj.originalText(chainLowerLimit:end),'Mass Value'),1,'first') + chainLowerLimit - 1;
                    massStr = obj.originalText{massInd};
                    startInd = strfind(massStr,'Actual="')+8;
                    endInd = strfind(massStr,'"/')-1;
                    mass = str2double(massStr(startInd:endInd))/1000;
                    obj.mass(i,1) = mass;

                    %keyboard

                    %Extract the overall dimensions for the segment
                    if i == 1
                        lengthInd = find(contains(obj.originalText(chainLowerLimit:end),'Length Value'),1,'last') + chainLowerLimit - 1;

                    else
                        lengthInd = find(contains(obj.originalText(chainLowerLimit:end),'Length Value'),obj.numBodies-(i-1),'first') + chainLowerLimit - 1;
                        lengthInd = lengthInd(end);
                    end

                    widthInd = lengthInd+1;
                    heightInd = widthInd+1;

                    lengthStr = obj.originalText{lengthInd};
                    widthStr = obj.originalText{widthInd};
                    heightStr = obj.originalText{heightInd};

                    startInd = strfind(lengthStr,'Actual="')+8;
                    endInd = strfind(lengthStr,'"/')-1;
                    l = str2double(lengthStr(startInd:endInd));
                    obj.lengths(i,1) = l;

                    startInd = strfind(widthStr,'Actual="')+8;
                    endInd = strfind(widthStr,'"/')-1;
                    width = str2double(widthStr(startInd:endInd));
                    obj.width(i,1) = width;

                    startInd = strfind(heightStr,'Actual="')+8;
                    endInd = strfind(heightStr,'"/')-1;
                    height = str2double(heightStr(startInd:endInd));
                    obj.height(i,1) = height;


                    %%% Now for the joint
                    if i > 1

                        %Move up the minimum line value now that that body is done.
                        jointFound = strfind(obj.originalText,['<Name>',obj.joints{i},'</Name>']);
                        nextJointInd = find(~cellfun(@isempty,jointFound),1,'first');
                        chainLowerLimit = nextJointInd;

                        %Find what kind of joint this is. If it is
                        %prismatic, it does not rotate anything. If it is a hinge, then we
                        %need to account for that rotation.
                        jointTypeInd = find(contains(obj.originalText(chainLowerLimit:end),'<PartType'),1,'first') + chainLowerLimit - 1;
                        jointTypeStr = obj.originalText(jointTypeInd);

                        %pull out the joint type
                        typeBegin = strfind(jointTypeStr,'.');
                        typeEnd = strfind(jointTypeStr,'<');
                        obj.jointTypes{i} = jointTypeStr{:}(typeBegin{:}(end)+1:typeEnd{:}(end)-1);

                        %Find the orientation of that joint
                        try jointRotInd = find(contains(obj.originalText(chainLowerLimit:end),'<Rotation'),1,'first') + chainLowerLimit - 1;
                        catch
                            error(['Joint ',obj.joints{i},' not found. Ensure that joints are listed correctly.']);
                        end
                        jointRotStr = obj.originalText(jointRotInd+1:jointRotInd+3);

                        for j=1:3
                            startInd = strfind(jointRotStr{j},'Actual="')+8;
                            endInd = strfind(jointRotStr{j},'"/>')-1;
                            obj.eulerAnglesJoints(j,i) = pi/180*str2double(jointRotStr{j}(startInd:endInd));
                        end

                        %Compute the rotation of the joint's axis within the local frame.
                        rotx_j = [1,0,0;0,cos(obj.eulerAnglesJoints(1,i)),-sin(obj.eulerAnglesJoints(1,i));0,sin(obj.eulerAnglesJoints(1,i)),cos(obj.eulerAnglesJoints(1,i))];
                        roty_j = [cos(obj.eulerAnglesJoints(2,i)),0,sin(obj.eulerAnglesJoints(2,i));0,1,0;-sin(obj.eulerAnglesJoints(2,i)),0,cos(obj.eulerAnglesJoints(2,i))];
                        rotz_j = [cos(obj.eulerAnglesJoints(3,i)),-sin(obj.eulerAnglesJoints(3,i)),0;sin(obj.eulerAnglesJoints(3,i)),cos(obj.eulerAnglesJoints(3,i)),0;0,0,1];

                        %Find the location of the joint within the frame
                        jointPosInd = find(contains(obj.originalText(chainLowerLimit:end),'<LocalPosition'),1,'first') + chainLowerLimit - 1;
                        jointPosStr = obj.originalText(jointPosInd+1:jointPosInd+3);

                        for j=1:3
                            startInd = strfind(jointPosStr{j},'Actual="')+8;
                            endInd = strfind(jointPosStr{j},'"/>')-1;
                            obj.posJoints(j,i) = str2double(jointPosStr{j}(startInd:endInd));
                        end

                        %Construct the initial exponential twist for the
                        %joint. Use the absolute rotation of the parent
                        %body and the local rotation of the joint to find
                        %its absolute rotation.

                        %Rotation of the joint's axis with respect to the
                        %parent's body.
                        R_local = rotx_j*roty_j*rotz_j;
                        R_parent = obj.gst0Bodies(1:3,1:3,i);

                        %Absolute rotation of the joint's axis.
                        R_abs = R_parent*R_local;

                        %Save this value into the exponential twist
                        %formula.
                        obj.gst0Joints(:,:,i) = eye(4);
                        obj.gst0Joints(1:3,1:3,i) = R_abs;

                        %Then, use the local rotation of the joint to
                        %transform the local position of the joint into the
                        %global frame, and add that to the parent body's
                        %location.
                        disp(obj.joints{i})
                        p_parent = obj.gst0Bodies(1:3,4,i);
                        q = obj.posJoints(:,i);
                        pJointAbs = p_parent + R_parent*q;

                        %Update the exponential twist formula.
                        obj.gst0Joints(1:3,4,i) = pJointAbs;

                        %Find the upper and lower rotation limits
                        jointLowerLimitInd = find(contains(obj.originalText(chainLowerLimit:end),'<LowerLimit'),1,'first') + chainLowerLimit - 1;
                        jointLowerLimitStr = obj.originalText{jointLowerLimitInd+3};

                        startInd = strfind(jointLowerLimitStr,'Actual="')+8;
                        endInd = strfind(jointLowerLimitStr,'"/>')-1;
                        obj.jointModelLowerLimits(i) = pi/180*str2double(jointLowerLimitStr(startInd:endInd));

                        jointUpperLimitInd = find(contains(obj.originalText(chainLowerLimit:end),'<UpperLimit'),1,'first') + chainLowerLimit - 1;
                        jointUpperLimitStr = obj.originalText{jointUpperLimitInd+3};

                        startInd = strfind(jointUpperLimitStr,'Actual="')+8;
                        endInd = strfind(jointUpperLimitStr,'"/>')-1;
                        obj.jointModelUpperLimits(i) = pi/180*str2double(jointUpperLimitStr(startInd:endInd));

                        %Check to see if a motor is enabled at this joint,
                        %making muscles unnecessary.
                        motorEnabledInd = find(contains(obj.originalText(chainLowerLimit:end),'<EnableMotor'),1,'first') + chainLowerLimit - 1;
                        motorEnabledStr = obj.originalText(motorEnabledInd);

                        motorEnabledCell = strfind(motorEnabledStr,'True');
                        %If motorEnabledCell is empty, then this joint
                        %cannot have a motor (ball and socket), and thus
                        %the motor cannot be enabled.
                        if isequal(size(motorEnabledCell),[1,1]) && ~isempty(motorEnabledCell{1}) && strcmp(obj.jointTypes{i},'Hinge')
                            obj.motorEnabled(i) = true;
                        else
                            %do nothing
                        end


                        %Move up the minimum line value now that that body is done.
                        try
                            nextJointInd = find(contains(obj.originalText,['<Name>',obj.bodies{i+1},'</Name>']),1,'first');
                            chainLowerLimit = nextJointInd;
                        catch
                            disp('End of chain reached.')
                        end
                    end
                end
                i = i + 1;
            end

            obj.jointHardwareLimits = [obj.jointModelLowerLimits,obj.jointModelUpperLimits];

            %% THIS SECTION NOT RELEVANT TO ROBOT

            if ~isRobot
                %%% Find muscles, their attachments, and associate them with
                %%% the proper bodies and joints.
                muscleFound = strfind(obj.originalText,'<Type>LinearHillMuscle</Type>');
                muscleInds = find(~cellfun(@isempty,muscleFound));
                muscleNameInds = muscleInds-2;

                %Does every joint in this chain have the same two character
                %prefix? For instance, LH or RF?
                onePrefix = true;
                for i=2:obj.numBodies-1
                    if ~strcmp(obj.joints{i}(1:2),obj.joints{i+1}(1:2))
                        onePrefix = false;
                    end
                end

                %If any two joints in this chain have different prefixes,
                %then query the user for the prefixes to use to categorize
                %this leg.
                if onePrefix
                    obj.legPrefix = obj.joints{2}(1:2);
                else
                    %leg_prefix = input('What is the prefix of objects belonging to this leg?\nA cell array of prefixes may be entered. ');
                    count = 0;
                    obj.legPrefix = nan;
                    while isnan(obj.legPrefix)
                        try
                            obj.legPrefix = obj.joints{end-count}(1:2);
                        catch
                            count = count+1;
                        end
                    end
                    warning('Not all joints in this leg use the same prefix. Continuing with the prefix of the last joint in this leg. If the last joint in your leg does not start with the correct Leg prefix, you will have problems')
                end
                obj.legPrefix

                %If there is only one prefix, then find the muscles whose
                %names begin with that prefix.
                if ischar(obj.legPrefix)
                    musclesForThisLeg = strfind(obj.originalText(muscleNameInds),['<Name>',obj.legPrefix]);
                elseif iscell(obj.legPrefix)
                    %If there are multiple prefixes, then check each prefix
                    %for muscles.
                    musclesForThisLeg = cell(length(muscleNameInds),1);
                    for i=1:length(obj.legPrefix)
                        tempMusclesForThisLeg = strfind(obj.originalText(muscleNameInds),['<Name>',obj.legPrefix{i}]);
                        for j=1:length(tempMusclesForThisLeg)
                            if ~isempty(tempMusclesForThisLeg{j})
                                musclesForThisLeg{j} = 1;
                            end
                        end
                    end
                end

                %Create a binary mask for the muscles in this organism that
                %belong to this leg.
                musclesForThisLegInds = cellfun(@isempty,musclesForThisLeg) == 0;

                %These are the indices of the names of muscles.
                muscleInds = muscleNameInds(musclesForThisLegInds);
                attachmentsToFind = cell(length(muscleInds),1);

                %Now that we know where the muscles are saved, we need to
                %extract all of the attachment IDs, find their names (to figure
                %out which body they belong to), and save their locations to
                %create joint objects.
                for i=1:length(muscleInds)
                    attachmentStartInd = find(contains(obj.originalText(muscleInds(i):end),'<Attachments>'),1,'first') + muscleInds(i) - 1;
                    attachmentEndInd = find(contains(obj.originalText(muscleInds(i):end),'</Attachments>'),1,'first') + muscleInds(i) - 1;

                    attachmentsToFind{i} = obj.originalText(attachmentStartInd + 1:attachmentEndInd - 1);
                    %remove the word "attach" from the strings
                    for j=1:length(attachmentsToFind{i})
                        attachmentsToFind{i}{j} = strrep(attachmentsToFind{i}{j},'Attach','');
                    end
                end

                attachmentNames = cell(size(attachmentsToFind));
                attachmentLocations = cell(size(attachmentsToFind));
                %find indices of the names of these attachments.
                for i=1:length(muscleInds)
                    for j=1:length(attachmentsToFind{i})
                        %save the names
                        idLoc = strfind(obj.originalText,attachmentsToFind{i}{j});
                        idInd = find(~cellfun(@isempty,idLoc),1,'first');
                        attachmentNames{i}{j} = idInd - 1;


                        %Find the position of that attachment. This mapping is
                        %identical to the names.
                        attachmentPositionInd = find(contains(obj.originalText(idInd:end),'<LocalPosition'),1,'first') + idInd - 1;
                        attachmentPositionStr = obj.originalText(attachmentPositionInd+1:attachmentPositionInd+3);

                        for k=1:3
                            startInd = strfind(attachmentPositionStr{k},'Actual="')+8;
                            endInd = strfind(attachmentPositionStr{k},'"/>')-1;
                            attachmentLocations{i}{j}(k,1) = str2double(attachmentPositionStr{k}(startInd:endInd));
                        end
                    end
                end

                %We will need to remove redundant attachments. If multiple muscles use the same
                %attachment, it will show up in our list twice. Therefore we
                %will start a list used_indices, and every time we add an
                %attachment's position to our record, we add the index of its
                %name to this list. Attachments whose indices already appear on
                %the list will be ignored.
                usedIndices = [];

                for i=2:obj.numBodies
                    obj.posAttachments{i} = cell(4,2);
                    %Find all instances of this joint's name within the
                    %attachment names, attach_names
                    for j=1:length(attachmentNames)
                        for k=1:length(attachmentNames{j})

                            %check to make sure the next attempt isn't already
                            %in this vector
                            if isempty(find(usedIndices == attachmentNames{j}{k},1,'first'))

                                jointNameFound = strfind(obj.originalText(attachmentNames{j}{k}),obj.joints{i});
                                %If we can identify which joint an attachment
                                %belongs to (i.e. anchors a muscle that acts across
                                %that joint), then save the position of that
                                %attachment in the appropriate location.
                                if ~cellfun(@isempty,jointNameFound) && isequal(size(attachmentLocations{j}{k}),[3,1])

                                    tempNameStr = char(obj.originalText(attachmentNames{j}{k}));
                                    tempNameStr = strrep(tempNameStr,'<Name>','');
                                    tempNameStr = strrep(tempNameStr,'</Name>','');

                                    %Determine if this attachment is associated
                                    %with an extensor or flexor by searching
                                    %for key phrases in its name.
                                    isFlexor = any([~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Flx')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'flx')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Flex')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'flex')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'flexor')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Flexor'))]);
                                    isExtensor = any([~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Ext')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'ext')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Extend')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'extend')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'extensor')),...
                                        ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Extensor'))]);
                                    keyboard
                                    %Determine if this attachment is on the
                                    %proximal or distal side of the joint it
                                    %actuates.
                                    %find bodies{i-1} or 'prox' in the name
                                    if iscell(obj.legPrefix)
                                        onProximalBool = zeros(length(obj.legPrefix)+4,1);
                                        for m=1:length(obj.legPrefix)
                                            onProximalBool(m) = ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),strrep(obj.bodies{i-1},strcat(obj.legPrefix{m},'_'),'')));
                                        end

                                        onProximalBool(end-3:end) = [~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Prox'));...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'prox'));...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Proximal'));...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'proximal'))];

                                        onProximal = any(onProximalBool);

                                    elseif ischar(obj.legPrefix)

                                        onProximal = any([~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),strrep(obj.bodies{i-1},strcat(obj.legPrefix,'_'),''))),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Prox')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'prox')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Proximal')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'proximal'))]);
                                    end

                                    %find bodies{i} or 'dist' in the name
                                    if iscell(obj.legPrefix)
                                        onDistalBool = zeros(length(obj.legPrefix)+6,1);
                                        for m=1:length(obj.legPrefix)
                                            onDistalBool(m) = ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),strrep(obj.bodies{i},strcat(obj.legPrefix{m},'_'),'')));
                                        end

                                        onDistalBool(end-5:end) = [~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Dis')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'dis')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Dist')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'dist')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Distal')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'distal'))];

                                        onDistal = any(onDistalBool);

                                    elseif ischar(obj.legPrefix)

                                        onDistal = any([~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),strrep(obj.bodies{i},strcat(obj.legPrefix,'_'),''))),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Dis')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'dis')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Dist')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'dist')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'Distal')),...
                                            ~cellfun(@isempty,strfind(obj.originalText(attachmentNames{j}{k}),'distal'))]);
                                    end
                                    %Save the position of the attachment in
                                    %question to the proper position in
                                    %obj.posAttachments. The first index is the
                                    %joint of interest. The second index is which
                                    %attachment it is, in this order: proximal
                                    %extensor, distal extensor, proximal flexor,
                                    %distal flexor.
                                    %Since attach_locs already holds position
                                    %vectors (in the local frame) of the
                                    %attachments indexed identically to the names,
                                    %we can map from the names (which we are
                                    %comparing above) to the positions (which we
                                    %are assigning below).
                                    if isExtensor
                                        if onProximal
                                            obj.posAttachments{i}{1,1} = [obj.posAttachments{i}{1,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{1,1},2)),') ');
                                            obj.posAttachments{i}{1,2} = [obj.posAttachments{i}{1,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        elseif onDistal
                                            obj.posAttachments{i}{2,1} = [obj.posAttachments{i}{2,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{2,1},2)),') ');
                                            obj.posAttachments{i}{2,2} = [obj.posAttachments{i}{2,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        else
                                            inputStr = strcat('The attachment:_ ',tempNameStr,' does not appear to belong to any body.\n',...
                                                'Which of the following does it belong to?\n');
                                            for m=1:obj.numBodies
                                                inputStr = strcat(inputStr,num2str(m),') ',bodies{m},'\n');
                                            end
                                            body = input(inputStr);
                                            if body == i - 1 %proximal
                                                obj.posAttachments{i}{1,1} = [obj.posAttachments{i}{1,1},attachmentLocations{j}{k}];
                                                numAttachments = strcat(num2str(size(obj.posAttachments{i}{1,1},2)),') ');
                                                obj.posAttachments{i}{1,2} = [obj.posAttachments{i}{1,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                            elseif body == i %distal
                                                obj.posAttachments{i}{2,1} = [obj.posAttachments{i}{2,1},attachmentLocations{j}{k}];
                                                numAttachments = strcat(num2str(size(obj.posAttachments{i}{2,1},2)),') ');
                                                obj.posAttachments{i}{2,2} = [obj.posAttachments{i}{2,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                            else
                                                disp('This class does not currently support multi-joint muscles.\nPlease select a body that directly interfaces with this joint.')
                                            end
                                        end
                                    elseif isFlexor
                                        if onProximal
                                            obj.posAttachments{i}{3,1} = [obj.posAttachments{i}{3,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{3,1},2)),') ');
                                            obj.posAttachments{i}{3,2} = [obj.posAttachments{i}{3,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        elseif onDistal
                                            obj.posAttachments{i}{4,1} = [obj.posAttachments{i}{4,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{4,1},2)),') ');
                                            obj.posAttachments{i}{4,2} = [obj.posAttachments{i}{4,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        else
                                            inputStr = strcat('The attachment:_ ',tempNameStr,' does not appear to belong to any body.\n',...
                                                'Which of the following does it belong to?\n');
                                            for m=1:obj.numBodies
                                                inputStr = strcat(inputStr,num2str(m),') ',bodies{m},'\n');
                                            end
                                            body = input(inputStr);
                                            if body == i - 1 %proximal
                                                obj.posAttachments{i}{3,1} = [obj.posAttachments{i}{3,1},attachmentLocations{j}{k}];
                                                numAttachments = strcat(num2str(size(obj.posAttachments{i}{3,1},2)),') ');
                                                obj.posAttachments{i}{3,2} = [obj.posAttachments{i}{3,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                            elseif body == i %distal
                                                obj.posAttachments{i}{4,1} = [obj.posAttachments{i}{4,1},attachmentLocations{j}{k}];
                                                numAttachments = strcat(num2str(size(obj.posAttachments{i}{4,1},2)),') ');
                                                obj.posAttachments{i}{4,2} = [obj.posAttachments{i}{4,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                            else
                                                disp('This class does not currently support multi-joint muscles.\nPlease select a body that directly interfaces with this joint.')
                                            end
                                        end

                                    else
                                        inputStr = ['The attachment:_',tempNameStr,' does not appear to belong to a flexor or extensor.\n',...
                                            'Which is it?\n1)Proximal extensor\n2)Distal extensor\n3)Proximal flexor\n4)Distal extensor\n5)Both flexor and extensor, proximal\n6)Both flexor and extensor, distal\n'];
                                        type = input(inputStr);
                                        if type == 1
                                            obj.posAttachments{i}{1,1} = [obj.posAttachments{i}{1,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{1,1},2)),') ');
                                            obj.posAttachments{i}{1,2} = [obj.posAttachments{i}{1,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        elseif type == 2
                                            obj.posAttachments{i}{2,1} = [obj.posAttachments{i}{2,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{2,1},2)),') ');
                                            obj.posAttachments{i}{2,2} = [obj.posAttachments{i}{2,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        elseif type == 3
                                            obj.posAttachments{i}{3,1} = [obj.posAttachments{i}{3,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{3,1},2)),') ');
                                            obj.posAttachments{i}{3,2} = [obj.posAttachments{i}{3,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        elseif type == 4
                                            obj.posAttachments{i}{4,1} = [obj.posAttachments{i}{4,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{4,1},2)),') ');
                                            obj.posAttachments{i}{4,2} = [obj.posAttachments{i}{4,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        elseif type == 5
                                            obj.posAttachments{i}{1,1} = [obj.posAttachments{i}{1,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{1,1},2)),') ');
                                            obj.posAttachments{i}{1,2} = [obj.posAttachments{i}{1,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];

                                            obj.posAttachments{i}{3,1} = [obj.posAttachments{i}{3,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{3,1},2)),') ');
                                            obj.posAttachments{i}{3,2} = [obj.posAttachments{i}{3,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        elseif type == 6
                                            obj.posAttachments{i}{2,1} = [obj.posAttachments{i}{2,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{2,1},2)),') ');
                                            obj.posAttachments{i}{2,2} = [obj.posAttachments{i}{2,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];

                                            obj.posAttachments{i}{4,1} = [obj.posAttachments{i}{4,1},attachmentLocations{j}{k}];
                                            numAttachments = strcat(num2str(size(obj.posAttachments{i}{4,1},2)),') ');
                                            obj.posAttachments{i}{4,2} = [obj.posAttachments{i}{4,2},numAttachments,obj.originalText{attachmentNames{j}{k}}(7:end-7),'\n'];
                                        else
                                            disp('This class does not currently support multi-muscle joints. Please select flexor or extensor.')
                                        end
                                    end
                                    %Flag attach_locs by replacing used values by -1.
                                    attachmentLocations{j}{k} = -1;

                                    %Add the index we just used to the exhausted
                                    %list.
                                    usedIndices = [usedIndices,attachmentNames{j}{k}]; %#ok<AGROW>
                                end
                            else
                                %We have already used this attachment
                            end
                        end
                    end
                end
                disp('Attachments located.')
                attachment_types = {' proximal extensor',' distal extensor',' proximal flexor',' distal flexor'};

                springFound = strfind(obj.originalText,'<Type>Spring</Type>');
                springInd = find(~cellfun(@isempty,springFound));
                springNameInd = springInd-2;
                springlRestInd = springInd+18;
                springkPassiveInd = springlRestInd + 1;
                springcPassiveInd = springkPassiveInd + 1;

                for i=1:length(springNameInd)
                    for j=2:obj.numBodies
                        if strfind(obj.originalText{springNameInd(i)},obj.joints{j})
                            springPresent(j) = true;
                            kPassivePosStr = obj.originalText(springkPassiveInd(i));
                            kPassiveStart = strfind(kPassivePosStr,'Actual="');
                            kPassiveStart = kPassiveStart{1}+8;
                            kPassiveEnd = strfind(kPassivePosStr,'"/>');
                            kPassiveEnd = kPassiveEnd{1} - 1;
                            kPassive(j)=str2double(kPassivePosStr{1}(kPassiveStart:kPassiveEnd));

                            cPassivePosStr = obj.originalText(springcPassiveInd(i));
                            cPassiveStart = strfind(cPassivePosStr,'Actual="');
                            cPassiveStart = cPassiveStart{1}+8;
                            cPassiveEnd = strfind(cPassivePosStr,'"/>');
                            cPassiveEnd = cPassiveEnd{1} - 1;
                            cPassive(j)=str2double(cPassivePosStr{1}(cPassiveStart:cPassiveEnd));
                            %  c_passive(j) = c_passive(j)*1000;
                            %1000 factor because Animatlab uses c_passive
                            %units in a dumb way. For some reason it is
                            %good now? Double check this always

                            lRestPosStr = obj.originalText(springlRestInd(i));
                            lRestStart = strfind(lRestPosStr,'Actual="');
                            lRestStart = lRestStart{1}+8;
                            lRestEnd = strfind(lRestPosStr,'"/>');
                            lRestEnd = lRestEnd{1} - 1;
                            lSpringRest(j)=str2double(lRestPosStr{1}(lRestStart:lRestEnd));
                        end
                    end
                end

                %Scott Edit
                %Find muscle properties. Compare to joints to find
                %what joint the muscle moves. Find whether it is
                %extensor or flexor. Find properties. Create
                %musprops for the joint
                muscleProps = NaN(2,8,length(joints));
                for i = 1:length(muscleInds)

                    %Find Kse
                    if i~=length(muscleInds)
                        ksePositionFound = strfind(obj.originalText(muscleInds(i):muscleInds(i+1)),'<Kse');
                        kpePositionFound = strfind(obj.originalText(muscleInds(i):muscleInds(i+1)),'<Kpe');
                        lWidthPosFound = strfind(obj.originalText(muscleInds(i):muscleInds(i+1)),'<Lwidth');
                        stimulusTensionPositionFound = strfind(obj.originalText(muscleInds(i):muscleInds(i+1)),'</StimulusTension>');
                    else
                        ksePositionFound = strfind(obj.originalText(muscleInds(i):end),'<Kse');
                        kpePositionFound = strfind(obj.originalText(muscleInds(i):end),'<Kpe');
                        lWidthPosFound = strfind(obj.originalText(muscleInds(i):end),'<Lwidth');
                        stimulusTensionPositionFound = strfind(obj.originalText(muscleInds(i):end),'</StimulusTension>');
                    end
                    ksePositionInd = find(~cellfun(@isempty,ksePositionFound))+muscleInds(i)-1;
                    ksePositionStr = obj.originalText(ksePositionInd);
                    kseStart = strfind(ksePositionStr,'Actual="');
                    kseStart = kseStart{1}+8;
                    kseEnd = strfind(ksePositionStr,'"/>');
                    kseEnd = kseEnd{1} - 1;
                    kse = str2double(ksePositionStr{1}(kseStart:kseEnd));

                    kpePositionInd = find(~cellfun(@isempty,kpePositionFound))+muscleInds(i)-1;
                    kpePositionStr = obj.originalText(kpePositionInd);
                    kpeStart = strfind(kpePositionStr,'Actual="');
                    kpeStart = kpeStart{1}+8;
                    kpeEnd = strfind(kpePositionStr,'"/>');
                    kpeEnd = kpeEnd{1} - 1;
                    kpe = str2double(kpePositionStr{1}(kpeStart:kpeEnd));

                    bPositionInd = kpePositionInd+1;
                    bPositionStr = obj.originalText(bPositionInd);
                    bStart = strfind(bPositionStr,'Actual="');
                    bStart = bStart{1}+8;
                    bEnd = strfind(bPositionStr,'"/>');
                    bEnd = bEnd{1} - 1;
                    b = str2double(bPositionStr{1}(bStart:bEnd));
                    %                     b=b/1000;
                    %                     %Animatlab treats damping weirdly... divide by 1000 to
                    %                     %put a bandaid on the boo-boo

                    lWidthPosInd = find(~cellfun(@isempty,lWidthPosFound))+muscleInds(i)-1;
                    lWidthPosStr = obj.originalText(lWidthPosInd);
                    lWidthStart = strfind(lWidthPosStr,'Actual="');
                    lWidthStart = lWidthStart{1}+8;
                    lWidthEnd = strfind(lWidthPosStr,'"/>');
                    lWidthEnd = lWidthEnd{1} - 1;
                    lWidth = str2double(lWidthPosStr{1}(lWidthStart:lWidthEnd));

                    stimulusTensionPositionInd = find(~cellfun(@isempty,stimulusTensionPositionFound))+muscleInds(i)-1;

                    amplitudePositionInd = stimulusTensionPositionInd - 3;
                    amplitudePositionStr = obj.originalText(amplitudePositionInd);
                    amplitudeStart = strfind(amplitudePositionStr,'Actual="');
                    amplitudeStart = amplitudeStart{1}+8;
                    amplitudeEnd = strfind(amplitudePositionStr,'"/>');
                    amplitudeEnd = amplitudeEnd{1} - 1;
                    amplitude = str2double(amplitudePositionStr{1}(amplitudeStart:amplitudeEnd));

                    steepnessPositionInd = stimulusTensionPositionInd - 2;
                    steepnessPositionStr = obj.originalText(steepnessPositionInd);
                    steepnessStart = strfind(steepnessPositionStr,'Actual="');
                    steepnessStart = steepnessStart{1}+8;
                    steepnessEnd = strfind(steepnessPositionStr,'"/>');
                    steepnessEnd = steepnessEnd{1} - 1;
                    steepness = str2double(steepnessPositionStr{1}(steepnessStart:steepnessEnd));
                    steepness = steepness/1000;
                    %Steepness is unitless, but animatlab multiplies
                    %steepness by xoffset-V in volts, whereas matlab we use
                    %mV. So, divide by 1000

                    xOffsetPositionInd = stimulusTensionPositionInd - 4;
                    xOffsetPositionStr = obj.originalText(xOffsetPositionInd);
                    xOffsetStart = strfind(xOffsetPositionStr,'Value="');
                    xOffsetStart = xOffsetStart{1}+7;
                    xOffsetEnd = strfind(xOffsetPositionStr,'" Scale');
                    xOffsetEnd = xOffsetEnd{1} - 1;
                    xOffset = str2double(xOffsetPositionStr{1}(xOffsetStart:xOffsetEnd));

                    yOffsetPositionInd = stimulusTensionPositionInd - 1;
                    yOffsetPositionStr = obj.originalText(yOffsetPositionInd);
                    yOffsetStart = strfind(yOffsetPositionStr,'Actual="');
                    yOffsetStart = yOffsetStart{1}+8;
                    yOffsetEnd = strfind(yOffsetPositionStr,'"/>');
                    yOffsetEnd = yOffsetEnd{1} - 1;
                    yOffset = str2double(yOffsetPositionStr{1}(yOffsetStart:yOffsetEnd));

                    %Determine which joint it is associated with
                    jointIndex = NaN(1,length(joints));
                    for j=1:length(joints)
                        if ~isempty(joints{j})
                            tempJointIndex = strfind(obj.originalText{muscleInds(i)},joints{j});
                            if ~isempty(tempJointIndex)
                                jointIndex(j) = tempJointIndex;
                            else
                                jointIndex(j) = 0;
                            end
                        else
                            jointIndex(j)=0;
                        end
                    end
                    jointIndex = find(jointIndex);

                    %Determine if this muscle is associated
                    %with an extensor or flexor by searching
                    %for key phrases in its name.
                    isFlexor = any([~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'Flx')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'flx')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'Flex')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'flex')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'flexor')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'Flexor'))]);
                    isExtensor = any([~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'Ext')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'ext')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'Extend')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'extend')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'extensor')),...
                        ~cellfun(@isempty,strfind(obj.originalText(muscleInds(i)),'Extensor'))]);

                    %Scott edit. Put it all into muscleProps
                    if isFlexor
                        muscleProps(1,:,jointIndex) = [kse, kpe, b, lWidth, amplitude, steepness, xOffset, yOffset];
                    elseif isExtensor
                        muscleProps(2,:,jointIndex) = [kse, kpe, b, lWidth, amplitude, steepness, xOffset, yOffset];
                    else
                        warning('This muscle is neither a flexor nor extensor. Properties are not being read correctly.')
                    end
                end

                %Scott edit. Find Neuron properties. Just need C, G, Er, Tonic stimulus.

                %Find indicies of all Neurons. cut out ones that we don't
                %want
                %Update 1/4/2017. Need more neurons now. Have Stiffness1,
                %Stiffness2, Intersection2, Integral, Integral2, Integral
                %Positive, Integral Negative, Zero.
                if ~isempty(obj.connectionMap)
                    allSynapses = [];
                    neuronFound = strfind(obj.originalText,'<ClassName>IntegrateFireGUI.DataObjects.Behavior.Neurons.NonSpiking</ClassName>');
                    neuronInds = find(~cellfun(@isempty,neuronFound));
                    neuronNameInds = neuronInds+27;
                    neuronNameStr = obj.originalText(neuronNameInds);
                    neuronsWeWant = strfind(neuronNameStr, obj.legPrefix);
                    for i=1:length(neuronsWeWant)
                        neuronsWeWant2(i) = ~isempty(neuronsWeWant{i}); %#ok<AGROW>
                    end
                    neuronNameStr = neuronNameStr(neuronsWeWant2);
                    neuronNameInds = neuronNameInds(neuronsWeWant2);


                    nprops = NaN(13,13,length(joints));
                    sprops = NaN(length(obj.connectionMap), 5, length(joints));

                    for i = 1:length(neuronNameInds)
                        %This finds all relevant neuron properties, one neuron
                        %at a time

                        if i~=length(neuronNameInds)
                            relativeSize_pos_found = strfind(obj.originalText(neuronNameInds(i):neuronNameInds(i+1)),'RelativeSize');
                        else
                            relativeSize_pos_found = strfind(obj.originalText(neuronNameInds(i):end),'RelativeSize');
                        end
                        relativeSize_pos_ind = find(~cellfun(@isempty,relativeSize_pos_found),1)+neuronNameInds(i)-1;
                        relativeSize_pos_str = obj.originalText(relativeSize_pos_ind);
                        relativeSize_start = strfind(relativeSize_pos_str,'Actual="');
                        relativeSize_start = relativeSize_start{1}+8;
                        relativeSize_end = strfind(relativeSize_pos_str,'"/>');
                        relativeSize_end = relativeSize_end{1} - 1;
                        RelativeSize=str2double(relativeSize_pos_str{1}(relativeSize_start:relativeSize_end));
                        G = RelativeSize; %in uS, nick szc says this is true, and he knows everything

                        Er_pos_ind = relativeSize_pos_ind-1;
                        Er_pos_str = obj.originalText(Er_pos_ind);
                        Er_start = strfind(Er_pos_str,'Value="');
                        Er_start = Er_start{1}+7;
                        Er_end = strfind(Er_pos_str,'" Scale');
                        Er_end = Er_end{1} - 1;
                        Er=str2double(Er_pos_str{1}(Er_start:Er_end));

                        timeConstant_pos_ind = relativeSize_pos_ind+1;
                        timeConstant_pos_str = obj.originalText(timeConstant_pos_ind);
                        timeConstant_start = strfind(timeConstant_pos_str,'Value="');
                        timeConstant_start = timeConstant_start{1}+7;
                        timeConstant_end = strfind(timeConstant_pos_str,'" Scale');
                        timeConstant_end = timeConstant_end{1} - 1;
                        timeConstant=str2double(timeConstant_pos_str{1}(timeConstant_start:timeConstant_end));

                        C = timeConstant*G;

                        tonic_stim_pos_ind = timeConstant_pos_ind + 22;
                        tonic_stim_pos_str = obj.originalText(tonic_stim_pos_ind);
                        tonic_stim_start = strfind(tonic_stim_pos_str,'Value="');
                        tonic_stim_start = tonic_stim_start{1}+7;
                        tonic_stim_end = strfind(tonic_stim_pos_str,'" Scale');
                        tonic_stim_end = tonic_stim_end{1} - 1;
                        tonic_stim=str2double(tonic_stim_pos_str{1}(tonic_stim_start:tonic_stim_end));

                        %Determine which joint it belongs to
                        jointIndex=NaN(1,length(joints));
                        for j = 1:length(joints)
                            if ~isempty(joints{j})
                                tempJointIndex = strfind(obj.originalText{neuronNameInds(i)},joints{j});
                                if ~isempty(tempJointIndex)
                                    jointIndex(j) = tempJointIndex;
                                else
                                    jointIndex(j) = 0;
                                end
                            else
                                jointIndex(j)=0;
                            end
                        end
                        jointIndex = find(jointIndex);

                        %Determine  which neurons this is by searching
                        %for key phrases in its name.
                        %Will need to update this everytime the network
                        %changes! :(
                        %Update 1/4/2017: added many neurons! yay
                        isFlexor = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Flx')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'flx')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Flex')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'flex')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'flexor')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Flexor'))]);
                        isExtensor = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Ext')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Extend')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'extend')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'extensor')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Extensor'))]);
                        is_PerceivedRotation = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Perceived')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Percieved'))]);
                        is_CommandedPosition = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Commanded')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Commanded Position'))]);
                        is_Intersection1 = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Intersection1')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'intersection1'))]);
                        is_Intersection2 = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Intersection2')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'intersection2'))]);
                        is_Stiffness1 = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Stiffness1')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'stiffness1'))]);
                        is_Stiffness2 = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Stiffness2')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'stiffness2'))]);
                        is_Integral1 = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Integral1')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'integral1'))]);
                        is_Integral2 = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Integral2')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'integral2'))]);
                        is_Zero = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Zero')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'zero'))]);
                        is_IntegralPostive = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Integral Positive')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Integral postive'))]);
                        is_IntegralNegative = any([~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Integral Negative')),...
                            ~cellfun(@isempty,strfind(obj.originalText(neuronNameInds(i)),'Integral negative'))]);

                        %Assigns node number to each neuron. This MUST match
                        %the structure of the network, and updated if more
                        %neurons are added
                        if isFlexor
                            neuron_node = 11;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif isExtensor
                            neuron_node = 12;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_PerceivedRotation
                            neuron_node = 13;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_CommandedPosition
                            neuron_node = 1;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_Intersection1
                            neuron_node = 2;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_Intersection2
                            neuron_node = 3;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_Stiffness1
                            neuron_node = 9;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_Stiffness2
                            neuron_node = 10;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_Integral1
                            neuron_node = 4;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_Integral2
                            neuron_node = 5;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_Zero
                            neuron_node = 6;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_IntegralPostive
                            neuron_node = 7;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        elseif is_IntegralNegative
                            neuron_node = 8;
                            nprops(neuron_node,:,jointIndex) = [C,G,NaN,NaN,Er,NaN,NaN,NaN,NaN,NaN,NaN,NaN,tonic_stim];
                        else
                            warning('neuron not identified. something is wrong, fix pls')
                        end

                        %While we are here, find the outgoing and ingoing synaptic connections.
                        if i~=length(neuronNameInds)
                            incomingSynapses_pos_found = strfind(obj.originalText(neuronNameInds(i):neuronNameInds(i+1)),'<InLinks>');
                        else
                            incomingSynapses_pos_found = strfind(obj.originalText(neuronNameInds(i):end),'<InLinks>');
                        end
                        incomingSynapses_pos_ind = find(~cellfun(@isempty,incomingSynapses_pos_found),1)+neuronNameInds(i)-1;

                        if i~=length(neuronNameInds)
                            outgoingSynapses_pos_found = strfind(obj.originalText(neuronNameInds(i):neuronNameInds(i+1)),'<OutLinks>');
                        else
                            outgoingSynapses_pos_found = strfind(obj.originalText(neuronNameInds(i):end),'<OutLinks>');
                        end
                        outgoingSynapses_pos_ind = find(~cellfun(@isempty,outgoingSynapses_pos_found),1)+neuronNameInds(i)-1;

                        if i~=length(neuronNameInds)
                            outgoingSynapses_End_pos_found = strfind(obj.originalText(neuronNameInds(i):neuronNameInds(i+1)),'</OutLinks>');
                        else
                            outgoingSynapses_End_pos_found = strfind(obj.originalText(neuronNameInds(i):end),'</OutLinks>');
                        end
                        outgoingSynapses_End_pos_ind = find(~cellfun(@isempty,outgoingSynapses_End_pos_found),1)+neuronNameInds(i)-1;

                        numberIncomingSynapses = outgoingSynapses_pos_ind-incomingSynapses_pos_ind-2;
                        numberOutgoingSynapses = outgoingSynapses_End_pos_ind - outgoingSynapses_pos_ind-1;
                        incomingSynapseID_all = [];
                        outgoingSynapseID_all = [];
                        if numberIncomingSynapses > 0
                            for k = 1:numberIncomingSynapses
                                incomingSynapseID_pos_found = incomingSynapses_pos_ind + k;
                                incomingSynapseID_pos_str = obj.originalText(incomingSynapseID_pos_found);
                                incomingSynapseID_start = strfind(incomingSynapseID_pos_str,'<ID>');
                                incomingSynapseID_start = incomingSynapseID_start{1}+4;
                                incomingSynapseID_end = strfind(incomingSynapseID_pos_str,'</ID>');
                                incomingSynapseID_end = incomingSynapseID_end{1} - 1;
                                incomingSynapseID=incomingSynapseID_pos_str{1}(incomingSynapseID_start:incomingSynapseID_end);
                                incomingSynapseID_all = [incomingSynapseID_all; incomingSynapseID];
                            end
                        end
                        if numberOutgoingSynapses > 0
                            for k = 1:numberOutgoingSynapses
                                outgoingSynapseID_pos_found = outgoingSynapses_pos_ind + k;
                                outgoingSynapseID_pos_str = obj.originalText(outgoingSynapseID_pos_found);
                                outgoingSynapseID_start = strfind(outgoingSynapseID_pos_str,'<ID>');
                                outgoingSynapseID_start = outgoingSynapseID_start{1}+4;
                                outgoingSynapseID_end = strfind(outgoingSynapseID_pos_str,'</ID>');
                                outgoingSynapseID_end = outgoingSynapseID_end{1} - 1;
                                outgoingSynapseID=outgoingSynapseID_pos_str{1}(outgoingSynapseID_start:outgoingSynapseID_end);
                                outgoingSynapseID_all = [outgoingSynapseID_all; outgoingSynapseID];
                            end
                        end
                        allSynapses = [allSynapses; incomingSynapseID_all; outgoingSynapseID_all];
                    end

                    %Scott Edit
                    %Here, we will read in Synapse Properties
                    %Connections are already found. First we will remove
                    %redundant connections. Then, we will remove adaptor
                    %connections, leaving us only with the synapses we are
                    %interested in. From there, we will figure out where it's
                    %connected and read in the properties

                    allSynapses = unique(allSynapses, 'rows');
                    m=1;
                    while m <= length(allSynapses(:,1))
                        connectionPosFound = strfind(obj.originalText(1:end),allSynapses(m,:));
                        connectionPosInd = find(~cellfun(@isempty,connectionPosFound));
                        for n = 1:length(connectionPosInd)
                            connection_name_pos_ind = connectionPosInd(n)-1;
                            isSynapse = any([~cellfun(@isempty,strfind(obj.originalText(connection_name_pos_ind),'synapse')),...
                                ~cellfun(@isempty,strfind(obj.originalText(connection_name_pos_ind),'Synapse'))]);
                            isAdapter = any([~cellfun(@isempty,strfind(obj.originalText(connection_name_pos_ind),'Adapter')),...
                                ~cellfun(@isempty,strfind(obj.originalText(connection_name_pos_ind),'adapter'))]);
                            if isAdapter
                                allSynapses(m,:) = [];
                                break
                            elseif isSynapse

                                %Save the position where the synapse ID can be used
                                %for later.

                                %Added "syn is normal" condition because in some
                                %cases, the synapses appears as OutLinks before
                                %InLinks
                                syn_is_normal = any([~cellfun(@isempty,strfind(obj.originalText(connectionPosInd(1)-2:connectionPosInd(1)),'<InLinks>')),...
                                    ~cellfun(@isempty,strfind(obj.originalText(connectionPosInd(1):connectionPosInd(1)+2),'</InLinks>'))]);
                                if any(syn_is_normal)
                                    synapseInPosInd(m) = connectionPosInd(1);
                                    synapseOutPosInd(m) = connectionPosInd(2);
                                else
                                    synapseInPosInd(m) = connectionPosInd(2);
                                    synapseOutPosInd(m) = connectionPosInd(1);
                                end
                                synapseTypePosInd(m) = connectionPosInd(3);
                                m=m+1;
                                break
                            end
                        end
                    end

                    %We now have the ID's of all of our synapses. For each
                    %synapse...
                    %1. Figure out it's pre and post synaptic neuron
                    %2. Figure out what joint it is associated with
                    %3. Figure out what synapse number it is in relation to connectionMap
                    %4. Read in it's properties
                    %5. Place into sprops
                    for m = 1:length(allSynapses(:,1))
                        %If large amounts of synapses are connected to one neuron, the 15 and
                        %25 may need to be changed. I assumed this range would
                        %be OK

                        %This will figure out the presynaptic Neuron
                        synapseOutLowerLimit = synapseOutPosInd(m)-25;
                        synapseOutUpperLimit = synapseOutPosInd(m)-15;
                        synapseOutPosFound = strfind(obj.originalText(synapseOutLowerLimit:synapseOutUpperLimit),'<Text>');
                        synapseOutPosIndCurrent = find(~cellfun(@isempty,synapseOutPosFound))+synapseOutLowerLimit-1;


                        syn_out_is_flexor = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Flx')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'flx')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Flex')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'flex')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'flexor')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Flexor'))]);
                        syn_out_is_extensor = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Ext')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Extend')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'extend')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'extensor')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Extensor'))]);
                        syn_out_is_PerceivedRotation = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Perceived')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Percieved'))]);
                        syn_out_is_CommandedPosition = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Commanded')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Commanded Position'))]);
                        syn_out_is_Intersection1 = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Intersection1')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'intersection'))]);
                        syn_out_is_Intersection2 = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Intersection2')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'intersection2'))]);
                        syn_out_is_Stiffness1 = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Stiffness1')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'stiffness'))]);
                        syn_out_is_Stiffness2 = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Stiffness2')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'stiffness2'))]);
                        syn_out_is_Integral1 = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Integral1')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'integral1'))]);
                        syn_out_is_Integral2 = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Integral2')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'integral2'))]);
                        syn_out_is_Zero = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Zero')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'zero'))]);
                        syn_out_is_IntegralPostive = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Integral Positive')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Integral postive'))]);
                        syn_out_is_IntegralNegative = any([~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Integral Negative')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapseOutPosIndCurrent),'Integral negative'))]);


                        %This will figure out the postsynaptic Neuron
                        synapse_in_lower_limit = synapseInPosInd(m)-25;
                        synapse_in_upper_limit = synapseInPosInd(m)-15;
                        synapse_in_pos_found = strfind(obj.originalText(synapse_in_lower_limit:synapse_in_upper_limit),'<Text>');
                        synapse_in_pos_ind_current = find(~cellfun(@isempty,synapse_in_pos_found))+synapse_in_lower_limit-1;

                        syn_in_is_flexor = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Flx')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'flx')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Flex')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'flex')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'flexor')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Flexor'))]);
                        syn_in_is_extensor = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Ext')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Extend')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'extend')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'extensor')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Extensor'))]);
                        syn_in_is_PerceivedRotation = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Perceived')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Percieved'))]);
                        syn_in_is_CommandedPosition = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Commanded')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Commanded Position'))]);
                        syn_in_is_Intersection1 = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Intersection1')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'intersection'))]);
                        syn_in_is_Intersection2 = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Intersection2')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'intersection2'))]);
                        syn_in_is_Stiffness1 = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Stiffness1')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'stiffness'))]);
                        syn_in_is_Stiffness2 = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Stiffness2')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'stiffness2'))]);
                        syn_in_is_Integral1 = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Integral1')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'integral1'))]);
                        syn_in_is_Integral2 = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Integral2')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'integral2'))]);
                        syn_in_is_Zero = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Zero')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'zero'))]);
                        syn_in_is_IntegralPostive = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Integral Positive')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Integral postive'))]);
                        syn_in_is_IntegralNegative = any([~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Integral Negative')),...
                            ~cellfun(@isempty,strfind(obj.originalText(synapse_in_pos_ind_current),'Integral negative'))]);

                        %Determine which joint it belongs to
                        jointIndex=NaN(1,length(joints));
                        for j = 1:length(joints)
                            if ~isempty(joints{j})
                                tempJointIndex = strfind(obj.originalText{synapse_in_pos_ind_current},joints{j});
                                if ~isempty(tempJointIndex)
                                    jointIndex(j) = tempJointIndex;
                                else
                                    jointIndex(j) = 0;
                                end
                            else
                                jointIndex(j)=0;
                            end
                        end
                        jointIndex = find(jointIndex);

                        %Find which synapse this is based on the pre and post
                        %synaptic neurons. This part will need to be edited if
                        %the ID's of the neurons change.


                        if syn_out_is_flexor
                            syn_out_id = 11;
                        elseif syn_out_is_extensor
                            syn_out_id = 12;
                        elseif syn_out_is_PerceivedRotation
                            syn_out_id = 13;
                        elseif syn_out_is_CommandedPosition
                            syn_out_id = 1;
                        elseif syn_out_is_Intersection1
                            syn_out_id = 2;
                        elseif syn_out_is_Intersection2
                            syn_out_id = 3;
                        elseif syn_out_is_Stiffness1
                            syn_out_id = 9;
                        elseif syn_out_is_Stiffness2
                            syn_out_id = 10;
                        elseif syn_out_is_Integral1
                            syn_out_id = 4;
                        elseif syn_out_is_Integral2
                            syn_out_id = 5;
                        elseif syn_out_is_Zero
                            syn_out_id = 6;
                        elseif syn_out_is_IntegralPostive
                            syn_out_id = 7;
                        elseif syn_out_is_IntegralNegative
                            syn_out_id = 8;
                        else
                            warning('synapse not properly not identified. something is wrong, fix pls')
                        end

                        if syn_in_is_flexor
                            syn_in_id = 11;
                        elseif syn_in_is_extensor
                            syn_in_id = 12;
                        elseif syn_in_is_PerceivedRotation
                            syn_in_id = 13;
                        elseif syn_in_is_CommandedPosition
                            syn_in_id = 1;
                        elseif syn_in_is_Intersection1
                            syn_in_id = 2;
                        elseif syn_in_is_Intersection2
                            syn_in_id = 3;
                        elseif syn_in_is_Stiffness1
                            syn_in_id = 9;
                        elseif syn_in_is_Stiffness2
                            syn_in_id = 10;
                        elseif syn_in_is_Integral1
                            syn_in_id = 4;
                        elseif syn_in_is_Integral2
                            syn_in_id = 5;
                        elseif syn_in_is_Zero
                            syn_in_id = 6;
                        elseif syn_in_is_IntegralPostive
                            syn_in_id = 7;
                        elseif syn_in_is_IntegralNegative
                            syn_in_id = 8;
                        else
                            warning('synapse not properly not identified. something is wrong, fix pls')
                        end

                        synapse_id_to_look_for = [syn_out_id syn_in_id];
                        [n_connections, ~] = size(obj.connectionMap);
                        for p = 1:n_connections
                            if obj.connectionMap(p,:) == synapse_id_to_look_for
                                synapse_id = p;
                                break
                            else
                                synapse_id = 0;
                            end
                        end

                        %Finally, figure out this synapses properties
                        syn_smooth_width = 0;

                        syn_type_ind = synapseTypePosInd(m)+39;
                        syn_type_id_str = obj.originalText(syn_type_ind);
                        syn_type_start = strfind(syn_type_id_str,'TypeID>');
                        syn_type_start = syn_type_start{1}+7;
                        syn_type_end = strfind(syn_type_id_str,'</Syn');
                        syn_type_end = syn_type_end{1} - 1;
                        syn_type=syn_type_id_str{1}(syn_type_start:syn_type_end);
                        syn_type = strcat('<ID>',syn_type);

                        syn_properties_found = strfind(obj.originalText(1:end),syn_type);
                        syn_properties_ind = find(~cellfun(@isempty,syn_properties_found),1);

                        E_syn_pos_found = strfind(obj.originalText(syn_properties_ind:syn_properties_ind+41),'<EquilibriumPotential');
                        E_syn_ind = find(~cellfun(@isempty,E_syn_pos_found))+syn_properties_ind-1;
                        G_syn_ind = E_syn_ind + 1;
                        E_lo_ind = G_syn_ind + 1;
                        E_hi_ind = E_lo_ind + 1;

                        E_syn_str = obj.originalText(E_syn_ind);
                        E_syn_start = strfind(E_syn_str,'Value="');
                        E_syn_start = E_syn_start{1}+7;
                        E_syn_end = strfind(E_syn_str,'" Scale');
                        E_syn_end = E_syn_end{1} - 1;
                        E_syn=str2double(E_syn_str{1}(E_syn_start:E_syn_end));

                        G_syn_str = obj.originalText(G_syn_ind);
                        G_syn_start = strfind(G_syn_str,'Actual="');
                        G_syn_start = G_syn_start{1}+8;
                        G_syn_end = strfind(G_syn_str,'"/>');
                        G_syn_end = G_syn_end{1} - 1;
                        G_syn=str2double(G_syn_str{1}(G_syn_start:G_syn_end));
                        G_syn = G_syn*10^6; %Put in uS

                        E_lo_str = obj.originalText(E_lo_ind);
                        E_lo_start = strfind(E_lo_str,'Value="');
                        E_lo_start = E_lo_start{1}+7;
                        E_lo_end = strfind(E_lo_str,'" Scale');
                        E_lo_end = E_lo_end{1} - 1;
                        E_lo=str2double(E_lo_str{1}(E_lo_start:E_lo_end));

                        E_hi_str = obj.originalText(E_hi_ind);
                        E_hi_start = strfind(E_hi_str,'Value="');
                        E_hi_start = E_hi_start{1}+7;
                        E_hi_end = strfind(E_hi_str,'" Scale');
                        E_hi_end = E_hi_end{1} - 1;
                        E_hi=str2double(E_hi_str{1}(E_hi_start:E_hi_end));


                        %Place into sprops
                        if synapse_id ~=0
                            sprops(synapse_id,:,jointIndex) = [G_syn,E_syn,E_lo,E_hi,syn_smooth_width];
                        else
                            warning('synapse not properly identified. fix.')
                        end

                    end
                else
                    %Do not read neuron/synapse properties

                end

            end
            %%
            %Now, we have all the information we need to construct a cell
            %of joint objects, which make up the leg.
            obj.jointObj = cell(obj.numBodies,1);

            for i=2:obj.numBodies

                if ~isRobot
                    %Make sure muscle attachments are listed in the proper
                    %order. If they are not, then our muscle force vector, as
                    %well as our muscle length, will be totally wrong.
                    %For each class of muscle attachment (proximal extensor,
                    %etc.)
                    for j=1:4
                        %If there is more than one attachment listed
                        if size(obj.posAttachments{i}{j,1},2) > 1
                            %Ask the user for them in order.
                            prompt_str = strcat('The',attachment_types{j},' for the_',obj.joints{i},' joint has multiple attachments.\n');
                            fprintf(prompt_str)
                            fprintf(obj.posAttachments{i}{j,2})
                            attach_order = input('Please input a vector arranging these names from the muscle terminus, \nalong the length of the muscle: ');
                            temp_mat = obj.posAttachments{i}{j,1};
                            obj.posAttachments{i}{j,1} = temp_mat(:,attach_order);
                        end
                    end
                else
                    obj.posAttachments{i}{1,1} = [];
                    obj.posAttachments{i}{2,1} = [];
                    obj.posAttachments{i}{3,1} = [];
                    obj.posAttachments{i}{4,1} = [];
                end

                %Instantiate the joint object.
                if isRobot
                    obj.jointObj{i} = joint(obj.posAttachments{i}{1,1},obj.posAttachments{i}{2,1},obj.posAttachments{i}{3,1},obj.posAttachments{i}{4,1},obj.posJoints(:,i),obj.eulerAnglesJoints(:,i),obj.posBodies(:,i),eulerAnglesBodies(:,i),(obj.ctrIndex == i),obj.motorEnabled(i),springPresent(i),[obj.jointModelLowerLimits(i),obj.jointModelUpperLimits(i)],obj.swingDuration,obj.maxStanceDuration,obj.maxDutyCycle,maxJointVelocity,joints{i},obj.isRobot,false,i,obj.legNumber);
                else
                    if isempty(obj.connectionMap)
                        if springPresent(i) == 1
                            obj.jointObj{i} = joint(obj.posAttachments{i}{1,1},obj.posAttachments{i}{2,1},obj.posAttachments{i}{3,1},obj.posAttachments{i}{4,1},obj.posJoints(:,i),obj.eulerAnglesJoints(:,i),obj.posBodies(:,i),eulerAnglesBodies(:,i),(obj.ctrIndex == i),obj.motorEnabled(i),springPresent(i),[obj.jointModelLowerLimits(i),obj.jointModelUpperLimits(i)],obj.swingDuration,obj.maxStanceDuration,obj.maxDutyCycle,maxJointVelocity,joints{i},obj.isRobot,false,i,obj.legNumber,muscleProps(:,:,i),lSpringRest(i),kPassive(i),cPassive(i));
                        else
                            obj.jointObj{i} = joint(obj.posAttachments{i}{1,1},obj.posAttachments{i}{2,1},obj.posAttachments{i}{3,1},obj.posAttachments{i}{4,1},obj.posJoints(:,i),obj.eulerAnglesJoints(:,i),obj.posBodies(:,i),eulerAnglesBodies(:,i),(obj.ctrIndex == i),obj.motorEnabled(i),springPresent(i),[obj.jointModelLowerLimits(i),obj.jointModelUpperLimits(i)],obj.swingDuration,obj.maxStanceDuration,obj.maxDutyCycle,maxJointVelocity,joints{i},obj.isRobot,false,i,obj.legNumber,muscleProps(:,:,i));
                        end
                    else
                        obj.jointObj{i} = joint(obj.posAttachments{i}{1,1},obj.posAttachments{i}{2,1},obj.posAttachments{i}{3,1},obj.posAttachments{i}{4,1},obj.posJoints(:,i),obj.eulerAnglesJoints(:,i),obj.posBodies(:,i),eulerAnglesBodies(:,i),(obj.ctrIndex == i),obj.motorEnabled(i),springPresent(i),[obj.jointModelLowerLimits(i),obj.jointModelUpperLimits(i)],obj.swingDuration,obj.maxStanceDuration,obj.maxDutyCycle,maxJointVelocity,joints{i},obj.isRobot,false,i,obj.legNumber,muscleProps(:,:,i),lSpringRest(i),kPassive(i),cPassive(i),nprops(:,:,i),loopParameters,sprops(:,:,i));
                    end
                end
            end
        end

        function [y,x1,xf,z1,zf] = legHeight(obj,config)
            tempToPlot = obj.toPlot;
            obj.setToPlot(false);
            [~,~,jointPts] = obj.computeJacobian(config,[]);
            y = jointPts(2,end);
            x1 = jointPts(1,1);
            xf = jointPts(1,end);
            xFTi = jointPts(1,end-2);
            zf = jointPts(3,end);
            z1 = jointPts(3,1);
            obj.setToPlot(tempToPlot);
        end

        function [J_base,footPoint,jointPoints,comPoints,omega,g_j,g_b,g_com,J_tool,figureHandle] = computeJacobian(obj,theta,figureHandle)
            if length(theta) == obj.numBodies - 2
                theta = [0;theta;0];
            else
                theta = zeros(obj.numBodies,1);
            end

            if isempty(figureHandle)
                tempToPlot = obj.toPlot;
                obj.setToPlot(false);
            end

            %For each joint, find its twist at the 0 configuration.
            if isempty(obj.twist)
                obj.twist = NaN(6,obj.numBodies);
                for i=2:obj.numBodies
                    %Use the initial position and orientation of each joint
                    %to compute q and omega in the 0 configuration. We can
                    %then transform these twists to compute the Jacobian.
                    q = obj.gst0Joints(1:3,4,i);
                    R = obj.gst0Joints(1:3,1:3,i);
                    o = R*[-1;0;0];
                    obj.twist(:,i) = [-cross(o,q);o];
                end
            end

            %Find the exponential twist matrix g for each frame/joint,
            %exp(xi*theta).
            g = NaN(4,4,obj.numBodies);
            for i=2:obj.numBodies
                g(:,:,i) = obj.g(obj.twist(:,i),theta(i));
            end

            %Find the location of each joint in space. Multiply the
            %cumulative exponential twist matrices up to that frame, and
            %postmultiply by the 0 configuration exponential twist to
            %obtain the rotation and translation of each joint.
            jointPoints = NaN(3,obj.numBodies);
            comPoints = NaN(3,obj.numBodies);
            g_j = NaN(4,4,obj.numBodies);
            g_b = NaN(4,4,obj.numBodies);
            g_com = NaN(4,4,obj.numBodies);
            g_cum = eye(4);
            for i=2:obj.numBodies
                g_cum = g_cum*g(:,:,i);
                g_j(:,:,i) = g_cum*obj.gst0Joints(:,:,i);
                g_b(:,:,i) = g_cum*obj.gst0Bodies(:,:,i);
                g_com(:,:,i) = g_cum*obj.gst0COMs(:,:,i);
                jointPoints(:,i) = g_j(1:3,4,i);
                comPoints(:,i) = g_com(1:3,4,i);
                bodyPoints(:,i) = g_b(1:3,4,i);
            end

            g_b(:,:,1) = obj.gst0Bodies(:,:,1);
            comPoints(:,1) = obj.gst0COMs(1:3,4,1);

            omega = NaN(3,obj.numBodies);
            for i=2:obj.numBodies
                omega(:,i) = g_j(1:3,1:3,i)*[-1;0;0];
            end

            %The jacobian matrix
            J_base = zeros(6,obj.numBodies-1);
            J_tool = zeros(6,obj.numBodies-1);

            for i=2:obj.numBodies
                if strcmp(obj.jointTypes{i},'Hinge')
                    J_tool(1:3,i-1) = -cross(omega(:,i),jointPoints(:,i));
                    J_tool(4:6,i-1) = omega(:,i);

                    J_base(1:3,i-1) = -cross(omega(:,i),jointPoints(:,i) - jointPoints(:,end));
                    J_base(4:6,i-1) = omega(:,i);
                elseif strcmp(obj.jointTypes{i},'Prismatic')
                    J_base(1:3,i-1) = omega(:,i);
                    J_tool(1:3,i-1) = omega(:,i);
                elseif strcmp(obj.jointTypes{i},'BallSocket')
                    %Ball-and-socket cannot be actuated, so its Jacobian
                    %might as well be 0.
                    %                     keyboard
                    %                     J_base(1:3,i-1) = omega(:,i);
                    %                     J_tool(1:3,i-1) = omega(:,i);

                else
                    error(['Joint type "',obj.jointTypes{i},'" unknown.'])
                end
            end
            %             %Because the tarsus joint is really used to find the point of
            %             %ground contact, we need to remove the column associated with
            %             %this joint.
            J_tool(:,1) = [];
            J_base(:,end) = [];
            jointPoints(:,1) = [];
            %             %             foot_point = joint_points(:,end) - joint_points(:,1);
            footPoint = jointPoints(:,end);
            omega(:,1) = [];

            if obj.toPlot

                figure(figureHandle)
                hold on
                plot3(jointPoints(1,:),-jointPoints(3,:),jointPoints(2,:),'Magenta','Linewidth',1)
                plot3(jointPoints(1,:),-jointPoints(3,:),jointPoints(2,:),'Magentao','Linewidth',1)
                plot3(comPoints(1,:),-comPoints(3,:),comPoints(2,:),'Blue*','Linewidth',1)
                plot3(bodyPoints(1,:),-bodyPoints(3,:),bodyPoints(2,:),'Redo','Linewidth',1)

                segmentLength = sqrt(sum(diff(jointPoints')'.^2));
                axisLen = .2*mean(segmentLength);
                obj.vecLen = mean(segmentLength);

                axisPoints = jointPoints + axisLen*omega;

                %Because we are using a joint at the end of the distal
                %segment for measuring the size of the last segment, we do
                %not need to draw its axis in the kinematic map.
                for i=1:obj.numBodies-2
                    x = [jointPoints(1,i),axisPoints(1,i)];
                    y = [jointPoints(2,i),axisPoints(2,i)];
                    z = [jointPoints(3,i),axisPoints(3,i)];
                    plot3(x,-z,y,'cyan')
                end

                grid on
                xlabel('x')
                ylabel('y')
                zlabel('z')
                view(0,0)
                axis equal
                drawnow
                hold off
            end

            if isempty(figureHandle)
                obj.setToPlot(tempToPlot);
            end

            %Remove the orientation of the last joint, which should not be
            %actuated.
            omega(:,end) = [];
        end

        function I = computeEquivalentInertia(obj,config,jointVec)
            %For each joint, we will find the location and orientation of
            %the distal masses. From this info we will find the moment of
            %inertia about the joint.

            if isempty(jointVec)
                jointVec = 2:obj.numBodies-1;
            else
                allJoints = 2:obj.numBodies-1;
                if isequal(union(allJoints,jointVec),allJoints)
                    %do nothing
                else
                    error(['Joints passed to Leg.computeEquivalentInertia must be within the set [2,',num2str(obj.numBodies-1),'].'])
                end
            end

            I = NaN(obj.numBodies-1,1);
            [~,~,~,~,g_j,g_b] = obj.computeJacobian(config,[]);

            %Find the COM of the segments distal of each joint i.
            for i=jointVec %for each joint
                p_j = g_j(1:3,4,i);
                R_j = g_j(1:3,1:3,i);
                R_x = obj.axis_angle_rotation([1;0;0],obj.eulerAnglesJoints(1,i));

                %Rotations of the joint about the x axis are not evident in
                %its configuration, because the x axis is the hinge axis.
                %Therefore we need to correct for that.
                R_j_inv = R_x*R_j';

                I_joint = zeros(3);
                for j=i:obj.numBodies-1 %for each body distal to this joint
                    %It's absolute rotation is known
                    R_b = g_b(1:3,1:3,j);

                    %Location of the body's center of mass, in the absolute
                    %frame.
                    p_b = g_b(1:3,4,j);

                    %vector pointing from the joint to the body's COM, in
                    %the joint's frame.
                    q = R_j_inv*(p_b - p_j);

                    %Moment of inertia of the COM as if it were a point
                    %mass (parallel axis thm). This is in the joint's
                    %frame.
                    I_com = obj.mass(j)*diag([q(2)^2 + q(3)^2,...
                        q(1)^2 + q(3)^2,...
                        q(1)^2 + q(2)^2]);

                    %We can find the relative rotation between the body's
                    %absolute frame, and the joint's frame. This is
                    %necessary to transform the moment of inertia from the
                    %body's frame (the principle axes) into the joint's
                    %frame.
                    R = R_j_inv*R_b;

                    %Moment of inertia of the RB about its COM in its own
                    %frame, transformed into the joint's frame of
                    %reference.
                    I_rb = R*obj.momInertia(:,:,j)*R';

                    %Inertia tensor for the segments distal of this joint,
                    %in the joint's frame.
                    I_joint = I_joint + I_com + I_rb;
                end

                %Use I_xx as the moment of inertia about this joint. This
                %is valid because the hinge joint counteracts torques that
                %occur in other directions, such as those produced by
                %off-diagonal inertia terms.
                obj.jointObj{i}.I_scalar = I_joint(1,1);
                I(i) = I_joint(1,1);

            end
            I(isnan(I)) = [];

        end

        function maxTorques = computeMaxTorques(obj,load)
            %We wish to find the max torque on each joint in the leg so we
            %can determine the maximum muscle tension, which is necessary
            %to tune the parameters of the muscle physical model.
            %Begin by initializing a vector to hold the maximum magnitude
            %of torque acting on the system.
            maxTorques = NaN(obj.numBodies-1,1);

            %only joints that can produce torque are useful in this
            %calculation.
            useableJoints = false(obj.numBodies,1);
            for i=2:obj.numBodies
                useableJoints(i) = obj.jointObj{i}.useable;
            end
            %             expected = [false;true(obj.numBodies-2,1);false];
            %
            %             if ~isequal(useable_joints,expected)
            %                 keyboard
            %                 error('When specifying the kinematic chain, the first joint must be empty, and the final joint must be the location of the foot. All other joints must have actuators, whether a motor or an antagonistic pair of muscles.')
            %             end


            %We are now going to find the maximum torque on each joint by a
            %force applied at the tarsus.

            %For each joint, we will find the maximum torque applied by a
            %unit force on the tarsus, directed perpendicular to the vector
            %from the hip to the foot.
            for i=2:obj.numBodies-1

                %Whether the user wants to plot results or not, we need to turn
                %it off temporarily so we can call computeJacobian thousands of
                %times in the optimizer without slowing/crashing the process.
                %Store it in a temporary variable, then change the property to
                %false.
                tempToPlot = get(obj,'toPlot');
                obj.setToPlot(false);

                %First, check if this joint is a hinge joint with muscles
                %acting on it.
                if obj.jointObj{obj.numBodies + 1 - i}.useable == 1
                    if i==2
                        %We can find the maximum torque trivially, by
                        %applying the maximum ground force perpendicularly
                        %to the last segment
                        [~,uForce] = obj.torqueWrapper([],1,i-1);
                        sol = 0;
                    else
                        %If we can work with this joint, define an objective
                        %that is the negative torque acting on the ith joint,
                        %squared. Minimizing this function will reveal the
                        %maximum torque.
                        f = @(x)obj.torqueWrapper(x,1,i-1);

                        %bnds changes every iteration, since our x vector grows
                        %with each. We are only controlling the i-2 most distal
                        %joints. Note that this means to find the maximum
                        %torque on the most distal joint, we do not need to
                        %change the robot configuration; rotating the force
                        %alone will reveal all possible torques.
                        bnds = obj.jointHardwareLimits(obj.numBodies+2-i:obj.numBodies-1,:);

                        x0 = bnds(:,2) - .01*range(bnds,2);

                        %Refine this solution with a bounded QN optimizer to
                        %make sure the point is optimal.
                        sol = fmincon(f,x0,[],[],[],[],bnds(:,1),bnds(:,2));

                        %Compute the unit vector of the force, and then scale
                        %it by the desired magnitude (probably some fraction of
                        %body weight).
                        [~,uForce] = obj.torqueWrapper(sol,1,i-1);
                    end
                    forceVec = load*uForce/2;

                    %Only now, when calculating the max torque, do we apply
                    %the value of the load. Up to this point, we are just
                    %using a unit force.
                    %Input the joint angles that we just computed, and
                    %insert them into the configuration of the leg.
                    currentConfig = zeros(obj.numBodies-2,1);
                    currentConfig(obj.numBodies + 1 - i:obj.numBodies-2) = sol;
                    [~,absFootPosition] = obj.computeJacobian(currentConfig,[]);
                    %Calculate J'*F to find the joint torques.
                    [J_base,~,jointPoints,jointAxes] = obj.computeJacobian(currentConfig,[]);
                    maxTorque = J_base'*[forceVec;0;0;0];

                    %We want to ensure that the resulting force is
                    %perpendiular to both the joint axis, and the vector
                    %from the joint of interest to the foot.
                    uJoint = jointAxes(:,obj.numBodies-i)/norm(jointAxes(:,obj.numBodies-i));
                    uLeg = (absFootPosition - jointPoints(:,obj.numBodies-i))/norm(absFootPosition - jointPoints(:,obj.numBodies-i));

                    if (uForce'*uLeg > sqrt(eps)) || (uForce'*uJoint > sqrt(eps))
                        error('Vectors are not perpendicular.')
                    end

                    %Reassign obj.toPlot to the user's preference.
                    obj.setToPlot(tempToPlot);

                    if obj.toPlot
                        kinematicFigure = figure;
                        obj.computeJacobian(currentConfig,kinematicFigure);

                        forceVecToDraw = [absFootPosition,absFootPosition + obj.vecLen*uForce];

                        figure(kinematicFigure)
                        hold on
                        plot3(forceVecToDraw(1,:),forceVecToDraw(2,:),forceVecToDraw(3,:),'black')
                        titleStr = strcat('Max torque configuration for joint ',num2str(obj.numBodies+1-i),'. Force vector: <',num2str(uForce(1)),',',num2str(uForce(2)),',',num2str(uForce(3)),'>');
                        title(titleStr)

                        drawnow
                        hold off
                    end

                    %store the maximum torque.
                    maxTorques(obj.numBodies + 1 - i) = abs(maxTorque(obj.numBodies - i));
                else
                    maxTorques(i) = NaN;
                end

                %Modify the joint object to save this max torque
                obj.jointObj{obj.numBodies + 1 - i}.setMaxTorque(maxTorques(obj.numBodies + 1 - i));
            end
        end

        function output = computeMaxMuscleForce(obj)
            %Now that we've calculated the maximum torque on each joint, we
            %can calculate the maximum force required in each muscle.
            applyToMuscles = 1;
            if applyToMuscles
                for i=2:obj.numBodies
                    if obj.jointObj{i}.useable
                        obj.jointObj{i}.computeMuscleForceWorstCase([]);
                    end
                end
            end
            output = 1;
        end

        function [torque_mag,force_vec] = torqueWrapper(obj,x,force,i)
            %Maximize the torque about a joint by manipulating the
            %orientation as well as the direction of the force.

            %The force vector has the magnitude of the 'load' passed by the
            %user to obj.computeMaxTorques. Rotate a unit vector by the
            %input angles.

            current_config = zeros(obj.numBodies-2,1);
            current_config(obj.numBodies-(length(x)+1):obj.numBodies-2) = x;
            [~,r_foot,rJoints,u_joints] = obj.computeJacobian(current_config,[]);

            r_joint = r_foot - rJoints(:,obj.numBodies-(1+i));
            u_segment = r_joint/norm(r_joint);
            u_joint = u_joints(:,obj.numBodies-(1+i));

            r_force = -cross(u_joint,u_segment);
            u_force = r_force/norm(r_force);

            force_vec = force*u_force;

            torque_mag = -norm(r_joint)*force;
        end

        function [synergy, kinematics, feasibleDistance, maxError] = linearFootTrajectorySynergy(obj,startPoint,movementVec,indsToModify,useJacobian)
            %Using a Moore-Penrose pseudo-inverse of the body-frame
            %Jacobian, we will find the joint velocities that move the foot
            %along desired trajectories with the least joint motion
            %possible.

            %Save the position of the root body for comparison to foot
            %positions
            bodyCenter = obj.posBodies(:,1);

            %If plotting, initialize a figure that will plot the desired
            %configurations of the leg.
            if obj.toPlot
                desiredLegConfigFig = figure('Position',[10,330,360,300]);
            end

            if any(isnan(movementVec))
                error('Movement vector must not include NaN.')
            end

            %If no initial foot position is provided, use the resting posture
            %saved in the leg object to start. Save the starting and ending
            %foot positions.
            if isempty(startPoint)
                startConfig = obj.getRestingConfig;
                if obj.toPlot
                    [~,~,rJoints] = obj.computeJacobian(startConfig,desiredLegConfigFig);
                else
                    [~,~,rJoints] = obj.computeJacobian(startConfig,[]);
                end
                %position relative to body's center
                startPoint = rJoints(:,end) - obj.posBodies(:,1);
                endPoint = startPoint + movementVec;

            elseif length(startPoint) == obj.numBodies-2
                startConfig = startPoint;
                if obj.toPlot
                    [~,~,rJoints] = obj.computeJacobian(startConfig,desiredLegConfigFig);
                else
                    [~,~,rJoints] = obj.computeJacobian(startConfig,[]);
                end
                %position relative to body's center
                startPoint = rJoints(:,end) - obj.posBodies(:,1);
                endPoint = startPoint + movementVec;
            else
                %do nothing
            end

            %We will find joint kinematics by starting from a given
            %configuration, and integrating J^-1*v_foot. We need to break
            %up the motion into small pieces, lest the process destabilize.
            %From preliminary tests, the foot velocity should not exceed
            %1mm/s. So, take 1000*s steps.

            s = norm(startPoint - endPoint);

            if useJacobian
                numSteps = ceil(s*1000);
                if numSteps <= 1
                    numSteps = 3;
                    s = numSteps/1000;
                    warning('Step length is being increased to 3 mm.')
                end
            else
                numSteps = ceil(s*50);
            end

            sVec = linspace(0,s,numSteps);
            kinematics = NaN(obj.numBodies-2,numSteps);
            deltaPosition = (endPoint - startPoint)/(numSteps-1);
            rFoot = NaN(3,numSteps);
            rFoot(:,1) = startPoint;
            positionError = NaN(1,numSteps);
            rCond = NaN(numSteps,1);
            rankVec = NaN(numSteps,1);

            kinematics = NaN(obj.numBodies-2,numSteps);
            rFoot = NaN(3,numSteps);
            rFoot(:,1) = startPoint;
            positionError = NaN(1,numSteps);
            rCond = NaN(numSteps,1);
            rankVec = NaN(numSteps,1);
            infeasibleFlag = false;

            if useJacobian
                %Compute the joint kinematics for the desired motion.
                finalStep = numSteps;
                i = 2;

                while i <= numSteps

                    prevConfig = kinematics(:,i-1);

                    if mod(i-2,10) == 0 && obj.toPlot
                        [J_base,~,rJoints] = obj.computeJacobian(prevConfig,desiredLegConfigFig);
                    else
                        toPlotTemp = obj.toPlot;
                        obj.setToPlot(false);
                        [J_base,~,rJoints] = obj.computeJacobian(prevConfig,[]);
                        obj.setToPlot(toPlotTemp);
                    end
                    rFoot(:,i) = rJoints(:,end) - bodyCenter;

                    inds = indsToModify-1;
                    allWithinBounds = false;
                    while ~allWithinBounds
                        J_pos = J_base(1:3,inds);
                        A = J_pos*J_pos';
                        rankVec(i) = rank(A);

                        rCond(i) = rcond(A);

                        b = deltaPosition;

                        nextConfig = prevConfig;

                        if rCond(i) < 1e-10 || rankVec(i) < 3
                            %The point is infeasible; we have too few degrees
                            %of freedom
                            allWithinBounds = false; %#ok<NASGU>
                            joints_outside_bounds = setdiff(indsToModify,inds+1);
                            error(['Desired kinematics for leg ',num2str(obj.legNumber),', joints ',num2str(joints_outside_bounds),', are infeasible. Starting kinematics should be changed.'])

                        else
                            configChange = J_pos'*linsolve(A,b);
                            nextConfig(inds) = nextConfig(inds) + configChange;

                            withinBounds = (nextConfig > obj.jointHardwareLimits(2:end-1,1)) &...
                                (nextConfig < obj.jointHardwareLimits(2:end-1,2));

                            allWithinBounds = all(withinBounds);
                        end

                        if ~allWithinBounds

                            %disp('out of bounds')
                            newInds = intersect(inds,find(withinBounds));
                            if isequal(newInds,inds) || isempty(inds)
                                %Point is infeasible
                                %disp('infeasible')
                                finalStep = i;
                                i = numSteps;
                                nextConfig = NaN(obj.numBodies-2,1);
                                allWithinBounds = true;
                            else
                                %disp('new inds')
                                inds = newInds;
                            end
                        end

                    end

                    kinematics(:,i) = nextConfig;

                    i = i + 1;
                end

                if obj.toPlot
                    obj.computeJacobian(nextConfig,desiredLegConfigFig);
                end
            else
                i = 2;
                toContinue = true;
                while i <= numSteps && toContinue
                    footPosition = startPoint + (i-1)*deltaPosition;
                    nearestConfig = obj.restingConfig(indsToModify-1);
                    f = @(x) obj.footPositionErrorDistance(x,footPosition,nearestConfig,indsToModify);
                    bnds = [obj.jointModelLowerLimits(indsToModify),...
                        obj.jointModelUpperLimits(indsToModify)];

                    x0 = kinematics(indsToModify-1,i-1);
                    [xf,~] = bnd_con_pen_bfgs(f,x0,bnds,[],[],[100,eps,1e-6],[],[],[],0);
                    f_final = f(xf);
                    positionError(i) = f_final(2);
                    if positionError(i) > 1e-3 %1 mm
                        toContinue = false;
                    end

                    kinematics(indsToModify-1,i) = xf;

                    if obj.toPlot
                        [~,~,rJoints] = obj.computeJacobian(xf,desiredLegConfigFig);
                    else
                        [~,~,rJoints] = obj.computeJacobian(xf,[]);
                    end
                    rFoot(:,i) = rJoints(:,end) - bodyCenter;

                    i = i + 1;
                end
                finalStep = i;
            end

            synergy = NaN(obj.numBodies-1,1);
            linearKinematics = NaN(obj.numBodies-2,finalStep-1);
            kinematics = kinematics(:,1:finalStep-1);
            rFoot = rFoot(:,1:finalStep-1);
            sVec = sVec(1:finalStep-1);
            feasibleDistance = sVec(end);

            if obj.toPlot
                fitLines = figure('Position',[370,690,360,300]);
                colors = {'blue','green','red','cyan','magenta','black'};
            end

            %for each joint, find a line of best fit. Start with a line
            %that has the slope of the endpoints, and an offset of the
            %initial position.
            for i=1:size(kinematics(:,1))
                [P,~] = polyfit(sVec,kinematics(i,:),1);
                m = P(1);
                b = P(2);

                outOfBoundsInds = union(find(kinematics(i,:) < obj.jointHardwareLimits(1+i,1)),...
                    find(kinematics(i,:) > obj.jointHardwareLimits(1+i,2)));

                if obj.toPlot
                    figure(fitLines)
                    hold on
                    grid on
                    plot(100*sVec,kinematics(i,:),colors{mod(i-1,6)+1},'Linewidth',2)
                    plot(100*sVec(outOfBoundsInds),kinematics(i,outOfBoundsInds),colors{mod(i-1,6)+1},'Linewidth',6)
                    plot(100*sVec,m*sVec + b,[colors{mod(i-1,6)+1},'--'],'Linewidth',3)
                    hold off
                end

                %The slope of the line gives us the proportional
                %contribution from each joint
                synergy(i) = m;

                %Now, find the foot positions of the linear "reconstructed"
                %kinematics for comparison.
                linearKinematics(i,:) = startConfig(i) + sVec*synergy(i);
            end

            if obj.toPlot
                figure(fitLines)
                hold on
                xlim([0,100*s])
                xlabel('Foot position (cm)')
                ylabel('Joint rotation (rad)')
                hold off

                act_fig = figure('Position',[370,330,360,300]);
            end

            %Initialize a matrix to store the reconstructed
            rFootLinear = NaN(3,finalStep-1);

            if useJacobian
                modNum = 10;
            else
                modNum = 1;
            end

            for i=1:finalStep-1
                if mod(i-1,modNum) == 0 && obj.toPlot
                    [~,~,rJoints] = obj.computeJacobian(linearKinematics(:,i),act_fig);
                else
                    toPlotTemp = obj.toPlot;
                    obj.setToPlot(false);
                    [~,~,rJoints] = obj.computeJacobian(linearKinematics(:,i),[]);
                    obj.setToPlot(toPlotTemp);
                end
                rFootLinear(:,i) = rJoints(:,end) - bodyCenter;
            end

            rFootDiff = rFoot - rFootLinear;
            positionError = sqrt(sum(rFootDiff.^2));
            maxError = max(positionError);

            if obj.toPlot
                figure('Position',[10,690,360,300]);
                hold on
                title('Euclidean foot position error')
                %                 plotyy(sVec,pos_error,sVec,norm_error)
                plot(100*sVec,100*positionError,'Linewidth',2)
                ylabel('cm')
                xlabel('Foot displacement from starting configuration')
                if ~isnan(maxError) && maxError ~= 0
                    axis([0,100*s,0,120*maxError])
                end
                grid on
                hold off

            end
        end

        function kinematics = prodExpInvKinematics(obj,startConfig,movementVec,numSteps,indsToModify,useJacobian)
            %If plotting, initialize a figure that will plot the desired
            %configurations of the leg.
            if obj.toPlot
                desiredConfigFigure = figure('Position',[10,330,360,300]);
            end

            if any(isnan(movementVec))
                error('Movement vector must not include NaN.')
            end
            %Normalize the movement vec so we can get all of the foot
            %points in the body frame
            movementVec = movementVec - movementVec(:,1);

            %Use the resting posture
            %saved in the leg object to start. Save the starting and ending
            %foot positions.
            if isempty(startConfig)
                startConfig = obj.getRestingConfig;
            end
            if obj.toPlot
                [~,~,rJoints,rCoMs] = obj.computeJacobian(startConfig,desiredConfigFigure);
                rTarsusCoM(:,1) = rCoMs(:,end-1);
            else
                [~,~,rJoints] = obj.computeJacobian(startConfig,[]);
            end
            %position relative to body's center
            startPoint = rJoints(:,end) - obj.posBodies(:,1);
            %Add this to the movementVec to find the foot point for all
            %steps
            rFoot = startPoint + movementVec;
            hold on
            if obj.toPlot
                plot3(rFoot(1,:),-rFoot(3,:),rFoot(2,:))
            end

            kinematics = NaN(obj.numBodies-2,numSteps);
            %Compute the joint kinematics for the desired motion.
            kinematics(:,1) = startConfig(:,1);
            if useJacobian
                %Compute the joint kinematics for the desired motion.
                finalStep = numSteps; %#ok<NASGU>
                i = 2;
                legNum = obj.legNumber;
                mobileJoints = find(obj.motorEnabled)-1;
                firstInd = mobileJoints(1);
                lastInd = mobileJoints(end);
                while i <= numSteps
                    prevConfig = kinematics(mobileJoints(1):mobileJoints(end),i-1);
                    if any(isnan(prevConfig))
                        keyboard
                    end
                    % opts = optimoptions('fmincon','StepTolerance',1e-10,'ConstraintTolerance',1e-10,'Display','final');
                    opts = optimoptions('fmincon','StepTolerance',1e-10,'ConstraintTolerance',1e-10,'Display','off');
                    % opts = optimoptions('fmincon','Display','final');
                    errorFun = @(vars) norm(vars - prevConfig);
                    nonLinCon = @(vars) footCond(obj,vars,rFoot(:,i),firstInd,lastInd);
                    if i >= 3
                        x0 = 2*prevConfig - kinematics(firstInd:lastInd,i-2);
                        [nextConfig,~,exitFlag,~] = fmincon(errorFun,x0,[],[],[],[],[],[],nonLinCon,opts); %
                        if exitFlag == -2
                            keyboard
                            [nextConfig,~,exitFlag,~] = fmincon(errorFun,2*prevConfig,[],[],[],[],[],[],nonLinCon,opts); %
                        end
                    else
                        [nextConfig,~,exitFlag,~] = fmincon(errorFun,2*prevConfig,[],[],[],[],[],[],nonLinCon,opts); %
                    end
                    kinematics(:,i) = zeros(obj.numBodies-2,1);
                    kinematics(mobileJoints(1):mobileJoints(end),i) = nextConfig;
                    if mod(i-2,3) == 0 && obj.toPlot
                        [J_base,~,rJoints,rCoMs] = obj.computeJacobian(kinematics(:,i),desiredConfigFigure);
                        keyboard
                    else
                        toPlotTemp = obj.toPlot;
                        obj.setToPlot(false);
                        [J_base,~,rJoints,rCoMs] = obj.computeJacobian(kinematics(:,i),[]);
                        obj.setToPlot(toPlotTemp);
                    end
                    rTarsusCoM(:,i) = rCoMs(:,end-1);
                    if exitFlag == -2
                        keyboard
                    end
                    i = i + 1;
                end
                mNew = (kinematics(:,5) - kinematics(:,end-4))/9;
                for p = 3:-1:0
                    kinematics(:,end-p) = kinematics(:,end-p-1) + mNew;
                end
                for p = 4:-1:1
                    kinematics(:,p) = kinematics(:,p+1) - mNew;
                end

            end
        end

        function [kinematics, infeasibleFlag] = moorePenroseInvKinematics(obj,startConfig,movementVec,numSteps,indsToModify,useJacobian)
            %Using a Moore-Penrose pseudo-inverse of the body-frame
            %Jacobian, we will find the joint velocities that move the foot
            %along desired trajectories with the least joint motion
            %possible.
            
            %Save the position of the root body for comparison to foot
            %positions
            bodyCenter = obj.posBodies(:,1);
            
            %If plotting, initialize a figure that will plot the desired
            %configurations of the leg.
            if obj.toPlot
                desiredConfigFigure = figure('Position',[10,330,360,300]);
            end
            
            if any(isnan(movementVec))
                error('Movement vector must not include NaN.')
            end
            
            %Use the resting posture
            %saved in the leg object to start. Save the starting and ending
            %foot positions.
            if isempty(startConfig)
                startConfig = obj.getRestingConfig;
                if obj.toPlot
                    [~,~,rJoints] = obj.computeJacobian(startConfig,desiredConfigFigure);
                else
                    [~,~,rJoints] = obj.computeJacobian(startConfig,[]);
                end
                %position relative to body's center
                startPoint = rJoints(:,end) - obj.posBodies(:,1);
            else
                if obj.toPlot
                    [~,~,rJoints] = obj.computeJacobian(startConfig,desiredConfigFigure);
                else
                    [~,~,rJoints] = obj.computeJacobian(startConfig,[]);
                end
                %position relative to body's center
                startPoint = rJoints(:,end) - obj.posBodies(:,1);
            end
            
            %We will find joint kinematics by starting from a given
            %configuration, and integrating J^-1*v_foot. We need to break
            %up the motion into small pieces, lest the process destabilize.
            %From preliminary tests, the foot velocity should not exceed
            %1mm/s. So, take 1000*s steps.
            
%             s = norm(startPoint - endPoint);
%             
%             if useJacobian
%                 numSteps = 200;
%             else
%                 numSteps = ceil(s*50);
%             end
            
            kinematics = NaN(obj.numBodies-2,numSteps);
            rFoot = NaN(3,numSteps);
            rFoot(:,1) = startPoint;
            positionError = NaN(1,numSteps);
            rCond = NaN(numSteps,1);
            rankVec = NaN(numSteps,1);
            infeasibleFlag = false;
            %Compute the joint kinematics for the desired motion.
            kinematics(:,1) = startConfig(:,1);
            if useJacobian
                %Compute the joint kinematics for the desired motion.
                finalStep = numSteps; %#ok<NASGU>
                i = 2;
                toContinue = true;
                
                while i <= numSteps && toContinue
                    
                    prevConfig = kinematics(:,i-1);
                    if any(isnan(prevConfig))
                        return;
                    end
                    
                    if mod(i-2,10) == 0 && obj.toPlot
                        [J_base,~,rJoints] = obj.computeJacobian(prevConfig,desiredConfigFigure);
                    else
                        toPlotTemp = obj.toPlot;
                        obj.setToPlot(false);
                        [J_base,~,rJoints] = obj.computeJacobian(prevConfig,[]);
                        obj.setToPlot(toPlotTemp);
                    end
                    rFoot(:,i) = rJoints(:,end) - bodyCenter;
                    
                    inds = indsToModify-1;
                    allWithinBounds = false;
                    while ~allWithinBounds
                        J_pos = J_base(1:3,inds);
                        A = J_pos*J_pos';
                        rankVec(i) = rank(A);
                        
                        rCond(i) = rcond(A);
                        
                        b = movementVec(:,i) - movementVec(:,i-1);
                        
                        nextConfig = prevConfig;

                        
                        if rCond(i) < 1e-3 || rankVec(i) < 3
                            keyboard
                            toContinue = false;

                            %The point is infeasible; we have too few degrees
                            %of freedom
%                             jointsOutsideBounds = setdiff(indsToModify,inds+1);
                            %                             error(['Desired kinematics for leg ',num2str(obj.legNumber),', joints ',num2str(joints_outside_bounds),', are infeasible. Starting kinematics should be changed.'])
                            infeasibleFlag = true;
                        else
                            configChange = J_pos'*linsolve(A,b);
                            try nextConfig(inds) = nextConfig(inds) + configChange;
                            catch
                                keyboard
                            end
                        end
                        
                        withinBounds = (nextConfig >= obj.jointHardwareLimits(2:end-1,1)) &...
                            (nextConfig <= obj.jointHardwareLimits(2:end-1,2));
                            
                        allWithinBounds = all(withinBounds);
                        
                        if ~allWithinBounds
                            
                            %disp('out of bounds')
                            try newInds = intersect(inds,find(withinBounds));
                            catch
                                keyboard
                            end
                            if isequal(newInds,inds) || isempty(inds)
                                %Point is infeasible
                                %disp('infeasible')
                                finalStep = i;
                                i = numSteps;
                                nextConfig = NaN(obj.numBodies-2,1);
                                allWithinBounds = true;
                            else
                                %disp('new inds')
                                inds = newInds;
                            end
                        end
                        
                    end
                    kinematics(:,i) = nextConfig;
                    i = i + 1;
                end
                
                if obj.toPlot
                    obj.computeJacobian(nextConfig,desiredConfigFigure);
                end
            else
                i = 2;
                toContinue = true;
                while i <= numSteps && toContinue
                    footPos = startPoint + (i-1)*deltaPosition;
                    nearestConfig = obj.restingConfig(indsToModify-1);
                    f = @(x) obj.footPositionErrorDistance(x,footPos,nearestConfig,indsToModify);
                    bnds = [obj.jointModelLowerLimits(indsToModify),...
                        obj.jointModelUpperLimits(indsToModify)];
                    
                    x0 = kinematics(indsToModify-1,i-1);
                    [xf,~] = bnd_con_pen_bfgs(f,x0,bnds,[],[],[100,eps,1e-6],[],[],[],0);
                    f_final = f(xf);
                    positionError(i) = f_final(2);
                    if positionError(i) > 1e-3 %1 mm
                        toContinue = false;
                    end
                    
                    kinematics(indsToModify-1,i) = xf;
                    
                    if obj.toPlot
                        [~,~,rJoints] = obj.computeJacobian(xf,desiredConfigFigure);
                    else
                        [~,~,rJoints] = obj.computeJacobian(xf,[]);
                    end
                    rFoot(:,i) = rJoints(:,end) - bodyCenter;
                    
                    i = i + 1;
                end
                finalStep = i;
            end
         
        end

        function [kinematicsResample,stKinematics,swKinematics,AEPID,PEPID,footPathResample] = steppingKinematics(obj,startConfig,jointVec,stanceVec,numSamps,stepHeight,stanceDuty,stepFrequency,stepProfile,terrainShape)
            if stanceDuty >= 1
                error('Stance phase duty must be < 1.')
            end

            colors = {'blue','green','red','cyan','magenta','black'};
            % tempToPlot = true;
            tempToPlot = false;

            useJacobian = true;
            if useJacobian
                modNum = 3;
            else
                modNum = 1;
            end

            obj.setToPlot(tempToPlot);

            %Generate the swing profile with 5th order polynomials
            q5 = @(t,p) p(1) + p(2)*t + p(3)*t.^2 + p(4)*t.^3 + p(5)*t.^4 + p(6)*t.^5;
            dq5 = @(t,p) p(2) + 2*p(3)*t + 3*p(4)*t.^2 + 4*p(5)*t.^3 + 5*p(6)*t.^4;
            ddq5 = @(t,p) 2*p(3) + 6*p(4)*t + 12*p(5)*t.^2 + 20*p(6)*t.^3;
            dddq5 = @(t,p) 6*p(4) + 24*p(5)*t + 60*p(6)*t.^2;
            
            % if obj.legNumber >= 5 %If it's a front limb, generate more in a teardrop shape
            %     stanceL = stanceVec(:,end)-stanceVec(:,1);
            %     posAdd = stanceL/3;
            %     stepNum = length(stanceVec);
            %     xMax = stanceVec(1,end);
            %     xMin = stanceVec(1,1);
            %     zMax = stanceVec(3,end);
            %     zMin = stanceVec(3,1);
            %     t0 = 0;
            %     tf = obj.swingDuration;
            %     tb1 = tf/2;
            %     tb2 = tf*3/4;
            % 
            %     A5_12 = [1 t0 t0^2 t0^3 t0^4 t0^5;...
            %         0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;...
            %         0 0 2 6*t0 12*t0^2 20*t0^3;...
            %         1 tb1 tb1^2 tb1^3 tb1^4 tb1^5;...
            %         0 1 2*tb1 3*tb1^2 4*tb1^3 5*tb1^4;...
            %         0 0 2 6*tb1 12*tb1^2 20*tb1^3];
            % 
            % 
            %     A5_23 = [1 tb1 tb1^2 tb1^3 tb1^4 tb1^5;...
            %         0 1 2*tb1 3*tb1^2 4*tb1^3 5*tb1^4;...
            %         0 0 2 6*tb1 12*tb1^2 20*tb1^3;...
            %         1 tb2 tb2^2 tb2^3 tb2^4 tb2^5;...
            %         0 1 2*tb2 3*tb2^2 4*tb2^3 5*tb2^4;...
            %         0 0 2 6*tb2 12*tb2^2 20*tb2^3];
            % 
            %     A5_34 = [1 tb2 tb2^2 tb2^3 tb2^4 tb2^5;...
            %         0 1 2*tb2 3*tb2^2 4*tb2^3 5*tb2^4;...
            %         0 0 2 6*tb2 12*tb2^2 20*tb2^3;...
            %         1 tf tf^2 tf^3 tf^4 tf^5;...
            %         0 1 2*tf 3*tf^2 4*tf^3 5*tf^4;...
            %         0 0 2 6*tf 12*tf^2 20*tf^3];
            % 
            %     A5_24 = [1 tb1 tb1^2 tb1^3 tb1^4 tb1^5;...
            %         0 1 2*tb1 3*tb1^2 4*tb1^3 5*tb1^4;...
            %         0 0 2 6*tb1 12*tb1^2 20*tb1^3;...
            %         1 tf tf^2 tf^3 tf^4 tf^5;...
            %         0 1 2*tf 3*tf^2 4*tf^3 5*tf^4;...
            %         0 0 2 6*tf 12*tf^2 20*tf^3];
            % 
            %     swingXVel = (xMin-xMax)/obj.swingDuration;
            %     swingXAcc = 0;
            % 
            %     swingZVel = (zMin-zMax)/obj.swingDuration;
            %     swingZAcc = 0;
            % 
            %     b5_12y = [0; 0; 0; stepHeight; 0; 0];
            %     b5_24y = [stepHeight; 0; 0; 0; 0; 0];
            % 
            %     b5_12x = [xMax; -swingXVel; swingXAcc; xMin; swingXVel; swingXAcc];
            %     b5_23x = [xMin; swingXVel; swingXAcc; xMin-posAdd(1); 0; 0];
            %     b5_34x = [xMin-posAdd(1); 0; 0; xMin; -swingXVel; swingXAcc];
            % 
            %     b5_12z = [zMax; -swingZVel; swingZAcc; zMin; swingZVel; swingZAcc];
            %     b5_23z = [zMin; swingZVel; swingZAcc; zMin - posAdd(3); 0; 0];
            %     b5_34z = [zMin - posAdd(3); 0; 0; zMin; -swingZVel; swingZAcc];
            % 
            %     t = linspace(t0,tf,stepNum);
            %     t12 = t(1:ceil(stepNum/2));
            %     t23 = t(ceil(stepNum/2)+1:ceil(stepNum*3/4));
            %     t34 = t(ceil(stepNum*3/4)+1:end);
            %     t24 = [t23 t34];
            % 
            %     p5_12y = linsolve(A5_12,b5_12y);
            %     p5_24y = linsolve(A5_24,b5_24y);
            % 
            %     p5_12x = linsolve(A5_12,b5_12x);
            %     p5_23x = linsolve(A5_23,b5_23x);
            %     p5_34x = linsolve(A5_34,b5_34x);
            % 
            %     p5_12z = linsolve(A5_12,b5_12z);
            %     p5_23z = linsolve(A5_23,b5_23z);
            %     p5_34z = linsolve(A5_34,b5_34z);
            % 
            %     y5_12 = q5(t12,p5_12y);
            %     y5_24 = q5(t24,p5_24y);
            %     y5_Sw = [y5_12 y5_24];
            % 
            %     x5_12 = q5(t12,p5_12x);
            %     x5_23 = q5(t23,p5_23x);
            %     x5_34 = q5(t34,p5_34x);
            %     x5_Sw = [x5_12 x5_23 x5_34];
            % 
            %     z5_12 = q5(t12,p5_12z);
            %     z5_23 = q5(t23,p5_23z);
            %     z5_34 = q5(t34,p5_34z);
            %     z5_Sw = [z5_12 z5_23 z5_34];
            % 
            %     %                 figure
            %     %                 plot3(x5_Sw,z5_Sw,y5_Sw,'-o')
            %     %                 keyboard
            % else %Otherwise us symmetrical profile

                stepNum = length(stanceVec);
                xMax = stanceVec(1,end);
                xMin = stanceVec(1,1);
                t0 = 0;
                tf = obj.swingDuration;
                tmid = tf/2;

                A5_13 = [1 t0 t0^2 t0^3 t0^4 t0^5;...
                    0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;...
                    0 0 2 6*t0 12*t0^2 20*t0^3;...
                    1 tf tf^2 tf^3 tf^4 tf^5;...
                    0 1 2*tf 3*tf^2 4*tf^3 5*tf^4;...
                    0 0 2 6*tf 12*tf^2 20*tf^3];
                A5_12 = [1 t0 t0^2 t0^3 t0^4 t0^5;...
                    0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;...
                    0 0 2 6*t0 12*t0^2 20*t0^3;...
                    1 tmid tmid^2 tmid^3 tmid^4 tmid^5;...
                    0 1 2*tmid 3*tmid^2 4*tmid^3 5*tmid^4;...
                    0 0 0 6 24*tmid 60*tmid^2];
                A5_23 = [1 tmid tmid^2 tmid^3 tmid^4 tmid^5;...
                    0 1 2*tmid 3*tmid^2 4*tmid^3 5*tmid^4;...
                    0 0 0 6 24*tmid 60*tmid^2;...
                    1 tf tf^2 tf^3 tf^4 tf^5;...
                    0 1 2*tf 3*tf^2 4*tf^3 5*tf^4;...
                    0 0 2 6*tf 12*tf^2 20*tf^3];

                swingXVel = (xMin-xMax)/obj.swingDuration;
                swingXAcc = 0;
                singleTimestep = (obj.swingDuration/(1-stanceDuty)*stanceDuty)/stepNum;
                swingYVel1 = (stanceVec(2,end)-stanceVec(2,end-1))/singleTimestep;
                swingYVel2 = (stanceVec(2,2)-stanceVec(2,1))/singleTimestep;

                b5_12y = [stanceVec(2,end); swingYVel1; 0; stepHeight; 0; 0];
                b5_23y = [stepHeight; 0; 0; stanceVec(2,1); swingYVel2; 0];
                b5_13x = [xMax; -swingXVel; swingXAcc; xMin; -swingXVel; swingXAcc];

                t = linspace(t0,tf,stepNum);
                t1 = t(1:(stepNum-1)/2);
                t2 = t((stepNum+1)/2:end);

                p5_1y = linsolve(A5_12,b5_12y);
                p5_2y = linsolve(A5_23,b5_23y);
                p5_x = linsolve(A5_13,b5_13x);

                y5_1 = q5(t1,p5_1y);
                y5_2 = q5(t2,p5_2y);
                y5_Sw = [y5_1 y5_2];
                x5_Sw = q5(t,p5_x);

                %                 Interpolate the z coordinates in swing based on stance

                if norm(x5_Sw) ~= 0
                    z5_Sw = interp1(stanceVec(1,:),stanceVec(3,:),x5_Sw,'pchip','extrap');
                else
                    z5_Sw = zeros(1,length(x5_Sw));
                end
                %                 figure
                %                 plot3(x5_Sw,z5_Sw,y5_Sw,'-o')
            % end

            footPathSwNew5th = [x5_Sw;y5_Sw;z5_Sw];

            footPathStFNew = stanceVec(:,(ceil(length(stanceVec)/2)):end);
            footPathStBNew = stanceVec(:,1:floor(length(stanceVec)/2));

            footPathNew5th = [footPathStFNew(:,1:end-1) footPathSwNew5th(:,1:end-1) footPathStBNew(:,1:end)];
            %Figure out what IDs AEP and PEP are for this foot trajectory
            if contains(terrainShape,'Flat')
                kinematicsNew5th = moorePenroseInvKinematics(obj,startConfig,footPathNew5th,length(footPathNew5th),jointVec,useJacobian);
            else
                kinematicsNew5th = prodExpInvKinematics(obj,startConfig,footPathNew5th,length(footPathNew5th),jointVec,useJacobian);
            end
            kinematics = kinematicsNew5th;

            if any(any(isnan(kinematics))) %If the Moore Penrose couldn't solve any of the configurations, it gives up and just leave NaN for the rest of them
                for k=1:length(kinematics)
                    if isnan(kinematics(1,k))
                        nanStartID = k;
                        break
                    end
                end
                %If this is the case, plot a visual of the footpath to show
                %where it failed, then error out
                figure
                plot3(footPathNew5th(1,1:(nanStartID-1)),footPathNew5th(3,1:(nanStartID-1)),footPathNew5th(2,1:(nanStartID-1)),'bo')
                hold on
                plot3(footPathNew5th(1,nanStartID),footPathNew5th(3,nanStartID),footPathNew5th(2,nanStartID),'rx')
                plot3(footPathNew5th(1,nanStartID+1:end),footPathNew5th(3,nanStartID+1:end),footPathNew5th(2,nanStartID+1:end),'k.')
                grid on
                legend('Solved','Error Point','Not Solved')
                error(['Full stepping kinematics could not be calculated.\n' ...
                    'Modify body translation or conditions in legHeightCon.m to ensure feasible motions.\n'...
                    'Solver failed at timestep: %s'], num2str(nanStartID))
            end

            if obj.toPlot
                hDesKin = figure;
                for i=1:size(kinematics,2)
                    if mod(i,modNum) == 0
                        [jacobianBodyFrame{i},~,jointPoints{i},comPoints{i},omegas{i},~,~] = obj.computeJacobian(kinematics(:,i),hDesKin); 
                    else
                        [jacobianBodyFrame{i},~,jointPoints{i},comPoints{i},omegas{i},~,~] = obj.computeJacobian(kinematics(:,i),[]); 
                    end
                    for j=1:6
                        jointPointsPlot{j}(:,i) = jointPoints{i}(:,j);
                    end
                end
            end
            %Scale the kinematics to the desired frequency and duty
            %Number of samples/steps in the data
            numSamps = size(kinematics,2);

            %Indices of the first stance and swing phases
            st1 = 1:length(footPathStFNew);
            sw1 = length(footPathStFNew)+1:length(footPathStFNew)+length(footPathSwNew5th);
            st2 = length(footPathStFNew)+length(footPathSwNew5th)+1:length(footPathNew5th);

            %We need to expand our kinematics to achieve the proper
            %stance/swing duty. Save the actual number of samples from
            %each.
            actualNumSwingSamps = length(sw1);
            actualNumStanceSamps = length(st1)+length(st2);

            %We will need to add samples to the stance phase to achieve the
            %proper proportion of stance and swing timesteps.
            additionalSamps = floor((stanceDuty*numSamps - actualNumStanceSamps)/(1-stanceDuty));
            %This assuming that you will never be moving fast enough to
            %need to remove stance phase points

            %We now have a desired number of samples, which will be used to
            %interpolate the kinematics.
            desiredNumStanceSamps = actualNumStanceSamps + additionalSamps;

            %Indices of the stance phase, now that we have resampled the
            %kinematics.
            st1Resample = linspace(st1(1),st1(end),floor(desiredNumStanceSamps/2));
            st2Resample = linspace(st2(1),st2(end),ceil(desiredNumStanceSamps/2));

            %Interpolate the stance phase kinematics
            st1Kinematics = NaN(size(kinematics,1),length(st1Resample));
            for i=1:size(kinematics,1)
                st1Kinematics(i,:) = interp1(st1,kinematics(i,st1),st1Resample);
            end
            for i=1:3
                footPathNew5thst1(i,:) = interp1(st1,footPathNew5th(i,st1),st1Resample); 
            end

            st2Kinematics = NaN(size(kinematics,1),length(st2Resample));
            for i=1:size(kinematics,1)
                st2Kinematics(i,:) = interp1(st2,kinematics(i,st2),st2Resample);
            end
            for i=1:3
                footPathNew5thst2(i,:) = interp1(st2,footPathNew5th(i,st2),st2Resample);
            end
            stKinematics = [st2Kinematics st1Kinematics];
            %Swing phase kinematics
            swKinematics = kinematics(:,length(footPathStFNew)+1:length(footPathStFNew)+length(footPathSwNew5th));
            footPathNew5thsw = footPathNew5th(:,length(footPathStFNew)+1:length(footPathStFNew)+length(footPathSwNew5th));
            %Full step, with stance and swing proportioned properly.
            kinematicsResample = [st1Kinematics,swKinematics,st2Kinematics];
            footPathResample = [footPathNew5thst1,footPathNew5thsw,footPathNew5thst2];
            trajDir = footPathResample(1,2)-footPathResample(1,1); %Pos if body moving backward
            if trajDir < 0
                [~,AEPID] = max(footPathResample(1,:)); %Here AEP means first point in stance
                [~,PEPID] = min(footPathResample(1,:));
            else
                [~,AEPID] = min(footPathResample(1,:)); %Here AEP means first point in stance
                [~,PEPID] = max(footPathResample(1,:));
            end
%             figure
%             for kinL= 1:length(kinematicsResample(:,1))
%                 hold on
%                 plot(kinematicsResample(kinL,:),'-o')
%             end
            %Redefine the num_samps to match the new number of samples.
            numSamps = size(kinematicsResample,2);
            period = 1/stepFrequency;

            %Make a time vector to produce the desired frequency.
            t = (1:numSamps)/numSamps*period;
            tStance = t(1:desiredNumStanceSamps);

            %Smooth out the disconnect between the start and end of the
            %kinematics due to error creep with a 3rd order polynomial
            %interpolation
            numPts = convergent(numSamps/6);
            if rem(numPts,2)
                numPts = numPts+1;
            end
            p1 = numPts/2 - 1;
            p2 = numPts/2;
            q3 = @(t,p) p(1) + p(2)*t + p(3)*t.^2 + p(4)*t.^3;
            t0 = 0;
            tf = numPts;
            tSmooth = linspace(t0,tf,numPts);

            A3_13 = [1 t0 t0^2 t0^3;...
                0 1 2*t0 3*t0^2;...
                1 tf tf^2 tf^3;...
                0 1 2*tf 3*tf^2];

            for k = 1:obj.numBodies-2
                y0 = kinematicsResample(k,end-p1);
                vy0 = kinematicsResample(k,end-p1) - kinematicsResample(k,end-(p1+1));
                yf = kinematicsResample(k,p2);
                vyf = kinematicsResample(k,p2+1)-kinematicsResample(k,p2);

                b3_13 = [y0; vy0; yf; vyf];
                p3 = linsolve(A3_13,b3_13);
                newKin = q3(tSmooth,p3);
                kinematicsResample(k,end-p1:end) = newKin(1:numPts/2);
                kinematicsResample(k,1:p2) = newKin(numPts/2+1:end);
            end

            %Combine three cycles into one vector.
            kinematicsComb = repmat(kinematicsResample,1,3);
            tComb = [t,period+t,2*period+t];
            tStanceComb = [tStance;period+tStance;2*period+tStance];

            kinematicsSmooth = NaN(size(kinematicsResample));
            kinematicsCombSmooth = NaN(size(kinematicsComb));

            for i=1:size(kinematicsComb,1)
                kinematicsSmooth(i,:) = smooth(kinematicsResample(i,:),50);
                kinematicsCombSmooth(i,:) = smooth(kinematicsComb(i,:),50);
            end

            if obj.toPlot
                hDesKin2 = figure;
                for i=1:size(kinematics,2)
                    if mod(i,modNum) == 0
                        obj.computeJacobian(kinematics(:,i),hDesKin2); %Just being used to plot here
                    end
                end

                figure;
                subplot(2,3,1:3);
                hold on
                for i=1:length(jointVec)
                    j = jointVec(i)-1;
                    try plot(tComb,kinematicsCombSmooth(j,:),colors{i})
                    catch
                        keyboard
                    end
                end

                %Plots joint kinematics for 3 steps to make sure everything
                %interpolates nicely
                plot(tStanceComb(1,:),zeros(size(tStanceComb(1,:)))+min(min(kinematicsResample))-.1*max(range(kinematicsResample)),'black','linewidth',2)
                plot(tStanceComb(2,:),zeros(size(tStanceComb(2,:)))+min(min(kinematicsResample))-.1*max(range(kinematicsResample)),'black','linewidth',2)
                plot(tStanceComb(3,:),zeros(size(tStanceComb(3,:)))+min(min(kinematicsResample))-.1*max(range(kinematicsResample)),'black','linewidth',2)
                grid on
                legend('1','2','3','4','stance')

                half_step_length = period*max(stanceDuty,1-stanceDuty);

                subplot(2,3,4)
                hold on
                for i=1:length(jointVec)
                    j = jointVec(i)-1;
                    plot(tComb(1:desiredNumStanceSamps),...
                        kinematicsCombSmooth(j,numSamps+(1:desiredNumStanceSamps)),colors{i})
                end
                title('Stance phase motions')
                legend('1','2','3','4')
                xlim([0,half_step_length])
                grid on

                subplot(2,3,5)
                hold on
                for i=1:length(jointVec)
                    j = jointVec(i)-1;
                    plot(tComb(1:actualNumSwingSamps),...
                        kinematicsCombSmooth(j,numSamps+desiredNumStanceSamps+(1:actualNumSwingSamps)),colors{i})
                end
                title('Swing phase motions')
                legend('1','2','3','4')
                xlim([0,half_step_length])
                grid on
            end
        end

        function success = assignProperties(obj)
            for i=2:obj.numBodies
                if obj.jointObj{i}.useable
                    obj.jointObj{i}.assign_props();
                end
            end

            if exist(obj.projectFile,'file') == 2
                success = 1;
                disp('New file created.')
            else
                success = 0;
            end
        end

        function output = trainControlNetwork(obj)
            for i=2:obj.numBodies
                if obj.jointObj{i}.useable
                    obj.jointObj{i}.trainPositionController();
                end
            end
            output = 1;
        end

        function angles = inverseCurrentToAngle(obj,currents)
            proportions = currents/20;
            angles = obj.jointHardwareLimits(:,1) + proportions.*(obj.jointHardwareLimits(:,2)-obj.jointHardwareLimits(:,1));
        end

        function success = computeRobotAdapterGains(obj,jointHardwareLimits,saveLoc)
            for i=2:obj.numBodies
                if obj.jointObj{i}.useable
                    obj.jointObj{i}.computeRobotAdapterGains(jointHardwareLimits,saveLoc);
                end
            end
            success = 1;
        end

        function success = setRestingConfig(obj,restingConfig,terrainShape)

            if length(restingConfig) == obj.numBodies-2
            else
                error('Resting configuration is not the proper length.')
            end
            if ~contains(terrainShape, 'Ball')
                for i = 1:obj.numBodies-2
                    if restingConfig(i) > obj.jointModelUpperLimits(i+1) || restingConfig(i) < obj.jointModelLowerLimits(i+1)
                        error('Resting Configuration of at least one joint is outside of the joint limits')
                    end
                end
            end

            obj.restingConfig = restingConfig;
            for i=1:obj.numBodies-2
                obj.jointObj{i+1}.setRestAngle(restingConfig(i));
            end

            success = true;
        end

        function config = getRestingConfig(obj)
            config = NaN(obj.numBodies-2,1);

            for i=2:obj.numBodies-1
                config(i-1) = obj.jointObj{i}.restingPostureAngle;
            end
        end

        function theta = angleBetweenSegments(obj,config,joint)
            numConfigs = size(config,2);
            theta = NaN(size(config,2),1);

            tempToPlot = obj.toPlot;
            obj.setToPlot(false);

            for i=1:numConfigs
                [~,~,pJoints] = obj.computeJacobian(config(:,i),[],false);
                rSegments = pJoints(:,2:end) - pJoints(:,1:end-1);

                u1 = -rSegments(:,joint-1)/norm(rSegments(:,joint-1));
                u2 = rSegments(:,joint)/norm(rSegments(:,joint));

                %asin returns [-pi/2,pi/2], so angles larger than pi/2 will
                %be returned as smaller than pi/2. Therefore, if the
                %vectors are pointed in opposite directions, we must use
                %the supplement of the asin angle.
                directionCosine = sign(u1'*u2);
                if directionCosine == 1
                    theta(i) = asin(norm(cross(u1,u2)));
                elseif directionCosine == -1
                    theta(i) = pi - asin(norm(cross(u1,u2)));
                end
            end

            obj.setToPlot(tempToPlot);
        end

        function toPlotVal = setToPlot(obj,toPlot)
            obj.toPlot = toPlot;
            for i=2:obj.numBodies
                obj.jointObj{i}.toPlot = toPlot;
            end
            toPlotVal = obj.toPlot;
        end

        function toPlot = getJointsToPlot(obj)
            toPlot = false(obj.numBodies-1,1);
            for i=2:obj.numBodies
                if obj.jointObj{i}.toPlot
                    toPlot(i-1) = true;
                end
            end
        end

        function propsCell = listAllLegParams(obj,legMirrored)

            propsCell = [];

            %Collect all of the low level controller properties
            for i=2:obj.numBodies
                if obj.jointObj{i}.useable
                    propsToConcat = obj.jointObj{i}.listAllJointParams(legMirrored);
                    propsCell = [propsCell;propsToConcat]; %#ok<AGROW>
                end
            end
        end

        function [cMatlab, cAnimatlab] = designAllDampersCockroach(obj,kVec)
            %Given a vector of the spring constant k for each joint, this
            %function will (hopefully!) calculate the spring damping
            %coefficient required for the joint to be at least critically
            %damped for all possible motion

            %Spring damping coefficient for a particular joint is ultimately a function  of the joint
            %angle theta of it's own joint and the joint angle of all
            %distal joints.

            %k_vec should be the same length of the joints vector, with NaN
            %on joints that are not actuated
            jointVec = find(~isnan(kVec(:,1)));
            xi = 1;
            cMatlab = NaN(length(kVec),1);
            cAnimatlab = NaN(length(kVec),1);
            anglesForMaxInertia = [];


            for i=flip(jointVec)'
                thetasCurrentJoint = linspace(obj.jointObj{i}.possible_rotations(1),obj.jointObj{i}.possible_rotations(end),50);
                cCurrent = zeros(50, 50);

                if i == jointVec(end)
                    %This is the most distal joint, and the c value is
                    %calculated a little differently
                    for j=1:length(thetasCurrentJoint)
                        thetaCurrent = thetasCurrentJoint(j);
                        [~,u] = obj.jointObj{i}.compute_muscle_length(thetaCurrent,'flx');
                        rRotated = axis_angle_rotation(obj.joint_obj{i},theta_current)*obj.joint_obj{i}.r_flx;
                        phi = norm(cross(rRotated,u))/(norm(rRotated)*norm(u));
                        rCurrent = norm(obj.jointObj{i}.r_flx);
                        kTorsionalCurrent = kVec(i)*rCurrent^2*sin(phi)^2;
                        configVector = [zeros((length(kVec)-3),1); thetaCurrent];
                        computeEquivalentInertia(obj,configVector,jointVec');
                        I_current = obj.jointObj{i}.I_scalar;
                        omega_n = sqrt(kTorsionalCurrent/I_current);
                        cTorsionalCurrent = xi*2*omega_n*I_current;
                        cCurrent(j) = cTorsionalCurrent/(rCurrent^2*sin(phi)^2);
                    end

                    cMatlab(i) = max(max(cCurrent));
                    cAnimatlab(i) = cMatlab(i)/1000;
                    thetasPreviousJoint = thetasCurrentJoint;
                else
                    %every other joint

                    for j = 1:length(thetasCurrentJoint)
                        for k=1:length(thetasPreviousJoint)
                            thetaCurrent = thetasCurrentJoint(j);
                            thetaPrevious = thetasPreviousJoint(k);
                            [~,u] = obj.jointObj{i}.compute_muscle_length(thetaCurrent,'flx');
                            rRotated = axis_angle_rotation(obj.joint_obj{i},theta_current)*obj.joint_obj{i}.r_flx;
                            phi = norm(cross(rRotated,u))/(norm(rRotated)*norm(u));
                            rCurrent = norm(obj.jointObj{i}.r_flx);
                            kTorsionalCurrent = kVec(i)*rCurrent^2*sin(phi)^2;

                            if all(jointVec(1:end-1)+1 == jointVec(2:end))
                                %Actuated joints are all consecutive. Add zeros to the config vector until we get a
                                %config vector of the desired length. This
                                %assumes the unactuated joints are all at the
                                %beginning of the kinematic chains (such as
                                %Body Joint)
                                configVector = [thetaCurrent;thetaPrevious;anglesForMaxInertia];
                                while length(configVector) ~= (length(kVec)-2)
                                    configVector = [0; configVector]; %#ok<AGROW>
                                end
                            else
                                configVector = [thetaCurrent;thetaPrevious;anglesForMaxInertia];
                                %first, insert zeros until we get a config
                                %vector that is one less that the desired
                                %length
                                while length(configVector) ~= (length(kVec)-3)
                                    configVector = [0; configVector]; %#ok<AGROW>
                                end
                                %Then, figure out where the unactuated
                                %joint is to add a zero
                                inBetween = find(jointVec(1:end-1)+1 ~= jointVec(2:end));
                                configVector = [configVector(1:inBetween);0;configVector(inBetween+1:end)];
                            end

                            obj.computeEquivalentInertia(configVector,jointVec');
                            I_current = obj.jointObj{i}.I_scalar;
                            omega_n = sqrt(kTorsionalCurrent/I_current);
                            cTorsionalCurrent = xi*2*omega_n*I_current;
                            cCurrent(j,k) = cTorsionalCurrent/(rCurrent^2*sin(phi)^2);
                        end
                    end
                    m = max(cCurrent);
                    [cMatlabNew, newAngleMaxInertiaInd] = max(m);
                    cMatlab(i) = cMatlabNew;
                    cAnimatlab(i) = cMatlab(i)/1000;

                    newAngleMaxInertia = thetasPreviousJoint(newAngleMaxInertiaInd);
                    anglesForMaxInertia = [anglesForMaxInertia;newAngleMaxInertia]; %#ok<AGROW>
                    thetasPreviousJoint = thetasCurrentJoint;
                end
            end
        end
    end %end methods

    methods(Static)
        function C = axis_angle_rotation(axis, angle)
            c = cos(angle);
            s = sin(angle);

            a1 = axis(1);
            a2 = axis(2);
            a3 = axis(3);

            C = [c+a1^2*(1-c), a1*a2*(1-c)-a3*s, a1*a3*(1-c)+a2*s;...
                a1*a2*(1-c)+a3*s, c+a2^2*(1-c), a2*a3*(1-c)-a1*s;...
                a1*a3*(1-c)-a2*s, a2*a3*(1-c)+a1*s, c+a3^2*(1-c)];

        end

        function g_mat = g( twist, angle )
            %G returns the transformation matrix for the given joint twist.

            g_mat = zeros(4);
            g_mat(4,4) = 1;

            if all(size(twist) == [6,1])
                %Provided twist is of the proper dimension
                vee = twist(1:3);
                omega = twist(4:6);

                if all(omega == 0)
                    %This is a sliding joint

                else
                    omegaSkew = leg.skew(omega);
                    expRot = eye(3) + sin(angle)*omegaSkew + (1-cos(angle))*omegaSkew^2;
                    g_mat(1:3,1:3) = expRot;
                    g_mat(1:3,4) = (eye(3)-expRot)*cross(omega, vee) + omega*omega'*vee*angle;
                end

            else
                error('Provided twist does not have the proper dimensions')
            end
        end

        function skewMat = skew( vec )
            %SKEW returns the skew matrix associated with a 3x1 vector.
            if all(size(vec) == [3,1])
                skewMat = [0 -vec(3), vec(2); vec(3), 0, -vec(1); -vec(2), vec(1), 0];
            else
                error('Vector is not the proper size to skew.')
            end
        end

        function ad_g = calculateAdjoint(g)
            %Given a rigid body transformation g, this function calculates the adjoint
            %of the mx

            R = g(1:3,1:3);
            p = g(1:3,4);
            pHat = skew(p);

            ad_g = [R, pHat*R; zeros(3), R];
        end

        function rbInv = rbInv(g)
            %Inverse of Rigid Body Transformation mx's
            R = g(1:3, 1:3);
            p = g(1:3, 4);
            rbInv = [R' -R'*p; 0 0 0 1];
        end
    end
end


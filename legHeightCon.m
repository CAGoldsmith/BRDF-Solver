function [c,ceq] = legHeightCon(theta,legObj,jointVec,desiredHeight,postureScaleFactor,terrainShape)

config = zeros(legObj.numBodies-2,1); %Angles of each joint in the leg (including immobile)
config(jointVec-1) = theta;
[height,x1,xf,z1,zf] = legObj.legHeight(config);
%x1, z1 are spatial coordinates of ThC joint
%xf, zf are spatial coordinates of tarsus tip
%height is current height of leg from origin
ceq(1) = height - desiredHeight; %Ensure that foot y pos creates desired height
c = [];
if legObj.legNumber <= 4 %If it's a middle or hind leg
    if legObj.legNumber == 1 %LH leg
        ceq(2) = (zf - z1) + 1.3*postureScaleFactor; 
        ceq(3) = (xf - x1) + 1.5*postureScaleFactor;
%         ceq(3) = (xf - x1) + 1.7*postureScaleFactor;
        if postureScaleFactor < 1 %Used only for Scale Drosophila to ensure fmincon doesn't break
            c(1) = config(1) - .2473;
        end
    elseif legObj.legNumber == 2 %RH leg
        ceq(2) = (zf - z1) - 1.3*postureScaleFactor;
        ceq(3) = (xf - x1) + 1.5*postureScaleFactor;
%         ceq(3) = (xf - x1) + 1.7*postureScaleFactor;
    elseif legObj.legNumber == 3 %LM leg
        if postureScaleFactor > .01
            if contains(terrainShape, 'Flat')
                ceq(2) = (zf - z1) + 2*postureScaleFactor;
            else
                ceq(2) = (zf - z1) + 2.5*postureScaleFactor;
            end
        end
        ceq(3) = (xf - x1) + .1*postureScaleFactor;
%         ceq(3) = (xf - x1);
    else %RM leg
        if postureScaleFactor > .01
            if contains(terrainShape, 'Flat')
                ceq(2) = (zf - z1) - 2*postureScaleFactor;
            else
                ceq(2) = (zf - z1) - 2.5*postureScaleFactor;
            end
        end
        ceq(3) = (xf - x1) + .1*postureScaleFactor;
%         ceq(3) = (xf - x1);
    end
else
    c(1) = (xf - x1) - 2.75*postureScaleFactor;
    if legObj.legNumber == 5 %LF
        ceq(2) = (zf - z1) + .65*postureScaleFactor;
    elseif legObj.legNumber == 6 %RF
        ceq(2) = (zf - z1) - .65*postureScaleFactor;
    end
end

end
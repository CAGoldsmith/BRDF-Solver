function CE = calculateCEeq(x,actuator,robot,omegas,n,TvarStart,TvarEnd,springK)
%Reconfigure the x vector to be in the form of the omega struct
for i=1:length(TvarStart)
    segNum(i) = robot.legObj{i}.numBodies-2;
    for j=1:segNum(i)
        if j ~= segNum(i)
            %omegas stored from most proximal to least proximal, while x is
            %created the opposite way. The struct will re-order torques to
            %follow the omega pattern
            xStruct{i}(:,j) = x(TvarEnd{i}-(3*(j-1)+2):TvarEnd{i}-(3*(j-1)));
        else
            xStruct{i}(:,j) = zeros(3,1);
        end
    end
end

CE = 0;

%Now cycle through each leg segment and add to the CE total if the joint is
%actuated
for i=1:length(TvarStart)
    for j=1:segNum(i)
        if robot.legObj{i}.motorEnabled(j+1)
            CE = CE + (dot(omegas{n}{i}(:,j),xStruct{i}(:,j)))^2;
        end
    end
end

CE = 1/(2*(actuator.k+springK))*CE;
end
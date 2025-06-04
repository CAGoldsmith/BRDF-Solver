        function [c,ceq] = kinematicsCondition(obj,thetas,rFoot)
            [~,~,rJoints] = obj.computeJacobian(thetas,[]);
            ceq = norm(rFoot - rJoints(:,end));
            c = [];
        end
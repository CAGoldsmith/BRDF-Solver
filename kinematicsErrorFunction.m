function error = kinematicsErrorFunction(obj,thetas,rFoot)
[~,~,rJoints] = obj.computeJacobian(thetas,[]);
error = norm(rFoot - rJoints(:,end));
end
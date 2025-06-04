function [c,ceq] = footCond(obj,vars,footPosDes,firstInd,lastInd)
c = [];
thetas = zeros(obj.numBodies-2,1);
thetas(firstInd:lastInd) = vars;
[~,~,rJoints] = obj.computeJacobian(thetas,[]);
footPosAct = rJoints(:,end);
ceq = footPosDes - footPosAct;
end
function error = errorCondProdExp(obj,vars,footPosDes,firstInd,lastInd)
thetas = zeros(obj.numBodies-2,1);
thetas(firstInd:lastInd) = vars;
[~,~,rJoints] = obj.computeJacobian(thetas,[]);
footPosAct = rJoints(:,end);
error = norm(footPosDes - footPosAct);
end
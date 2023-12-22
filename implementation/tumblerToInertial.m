function outputRef = tumblerToInertial(outputRef, outputEst)
% DESCRIPTION: Given the current estimated target orientation, and
% reference in tumbler frame obtain reference in LVLH frame
%
% INPUT:
%    outputRef           trajectory in tumbler frame from reference generator
%    outputEst           current estimated target orientation and parameters
%
% OUTPUT:
%	 outputRef  	     trajectory with additional data in LVLH frame
    tt = outputRef.time;
    zRefT = outputRef.xRefT;
    vRefT = outputRef.vRefT;
    zRefI = zeros(size(zRefT));
    vRefI = zeros(size(vRefT));
    if isscalar(tt)
        q = outputEst.tumblerState(4:7);
        zRefI(1:3,1) = quat.rotate(zRefT(1:3,1),q);
        zRefI(4:6,1) = quat.rotate(zRefT(4:6,1),q);
        vRefI(:,1) = quat.rotate(vRefT(:,1),q);
        outputRef.xRefI = zRefI;
        outputRef.vRefI = vRefI;
        return;
    end
    

    x0 = outputEst.tumblerState;

    param_mod = outputEst.parameters;
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~,stateRotTrue] = ode45(@(t,x) dynamics.odetumbler(t,x,param_mod), ...
                    tt,x0, options);
    stateRotTrue = stateRotTrue';

    for i = 1:length(tt)
        q = stateRotTrue(4:7,i);
        zRefI(1:3,i) = quat.rotate(zRefT(1:3,i),q);
        zRefI(4:6,i) = quat.rotate(zRefT(4:6,i),q);
        vRefI(:,i) = quat.rotate(vRefT(:,i),q);
    end

    outputRef.xRefI = zRefI;
    outputRef.vRefI = vRefI;

end


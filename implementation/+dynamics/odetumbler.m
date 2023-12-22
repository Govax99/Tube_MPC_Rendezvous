function dx = odetumbler(time, state, parameters)
%DESCRIPTION state-space form of tumbling-target dynamics for ode integration
% (linear dynamics in L-frame, angular dynamics in T-frame)
%
% INPUT:
%    time                current integration time
%    state               current chaser state [q_LT, w_IT_T]
%    parameters          structure containing parameters for the dynamics
%
% OUTPUT:
%	 dx  	             derivative of state
%
% LEGEND FRAMES:
%    L-frame             LVLH, local vertical-local horizontal frame, center on target cm
%    T-frame             relative frame, attached to target object, center on target cm
%
% MONOGRAM NOTATION:
%    p_XY_Z              p: symbol of physical quantity,
%                        X: "measured from",
%                        Y: target point/frame,
%                        Z: "expressed in" frame
[~,J_T,~,OM] = dynamics.set_parameters(parameters);
OM_IL_L = [0; 0; OM];

w_IT_T = state(1:3);
q_LT = state(4:7)./norm(state(4:7));

R_LT = quat.quat2rotm(q_LT);
w_LT_L = R_LT*w_IT_T - OM_IL_L;

A_T = quat.quat_kin_matrix(w_LT_L);

dx = [
    J_T\(-cross(w_IT_T, J_T*w_IT_T)); ...
    1/2*A_T*q_LT; ...
    ];
end

